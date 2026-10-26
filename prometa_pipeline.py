"""Download, QC, trim and quantify SRA runs listed in metadata_geral.csv.

Expected metadata columns (the existing table format): sample_id, condition, study_id.
Each row represents one SRA run. Library layout is detected from fasterq-dump output.
"""

from __future__ import annotations

import logging
import os
import shutil
import subprocess
from pathlib import Path

import pandas as pd


# Edit paths and resource settings here. PROMETA_METADATA can point to the existing CSV
# without changing its columns (for example, the CSV can remain in the WSL project).
ROOT = Path(__file__).resolve().parent
METADATA = Path(os.environ.get("PROMETA_METADATA", str(ROOT / "metadata_geral.csv"))).expanduser()
REFERENCE_INDEX = ROOT / "reference" / "salmon_index"
THREADS = 6
KEEP_INTERMEDIATES = True

SRA_DIR = ROOT / "sra_bruto"
RAW_DIR = ROOT / "raw_data"
TRIM_DIR = RAW_DIR / "trimmed"
QUANT_DIR = ROOT / "quant"
LOG_DIR = ROOT / "logs"
QC_DIR = ROOT / "qc"
TMP_DIR = ROOT / "tmp"
ENRICHED_METADATA = ROOT / "metadata_processado.csv"
RUN_STATUS = ROOT / "run_status.csv"


def _setup_logging() -> None:
    LOG_DIR.mkdir(parents=True, exist_ok=True)
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s [%(levelname)s] %(message)s",
        handlers=[logging.FileHandler(LOG_DIR / "pipeline.log", encoding="utf-8"), logging.StreamHandler()],
        force=True,
    )


def _load_metadata() -> pd.DataFrame:
    if not METADATA.is_file():
        raise FileNotFoundError(
            f"Não encontrei {METADATA}. Copie o metadata_geral.csv para a pasta do projeto "
            "ou altere METADATA no início de prometa_pipeline.py."
        )
    meta = pd.read_csv(METADATA, dtype=str).fillna("")
    required = {"sample_id", "condition", "study_id"}
    missing = required - set(meta.columns)
    if missing:
        raise ValueError(f"Colunas obrigatórias ausentes: {sorted(missing)}")
    meta["sample_id"] = meta["sample_id"].str.strip()
    meta["condition"] = meta["condition"].str.strip().str.lower()
    meta["study_id"] = meta["study_id"].str.strip()
    if (meta[["sample_id", "condition", "study_id"]] == "").any().any():
        raise ValueError("Há campos vazios em sample_id, condition ou study_id.")
    if meta["sample_id"].duplicated().any():
        duplicates = meta.loc[meta["sample_id"].duplicated(False), "sample_id"].tolist()
        raise ValueError(f"sample_id duplicado: {duplicates}")
    if not meta["condition"].isin(["control", "treatment"]).all():
        bad = sorted(set(meta.loc[~meta["condition"].isin(["control", "treatment"]), "condition"]))
        raise ValueError(f"condition deve ser control/treatment; valores encontrados: {bad}")
    for study, group in meta.groupby("study_id"):
        counts = group["condition"].value_counts()
        if set(counts.index) != {"control", "treatment"}:
            logging.warning("Estudo %s não contém as duas condições: %s", study, counts.to_dict())
        elif min(counts) < 2:
            logging.warning("Estudo %s tem menos de 2 runs em uma condição; DESeq2 não será confiável.", study)
    return meta


def _require_tools() -> None:
    required = ["prefetch", "fasterq-dump", "trim_galore", "fastqc", "salmon"]
    missing = [tool for tool in required if shutil.which(tool) is None]
    if missing:
        raise EnvironmentError(f"Ferramentas não encontradas no PATH: {', '.join(missing)}")
    if not REFERENCE_INDEX.is_dir():
        raise FileNotFoundError(f"Índice do Salmon não encontrado: {REFERENCE_INDEX}")


def _run(cmd: list[str], log_path: Path, cwd: Path = ROOT) -> None:
    log_path.parent.mkdir(parents=True, exist_ok=True)
    with log_path.open("a", encoding="utf-8") as log:
        log.write("\n$ " + subprocess.list2cmdline(cmd) + "\n")
        log.flush()
        completed = subprocess.run(cmd, cwd=cwd, stdout=log, stderr=subprocess.STDOUT, text=True)
    if completed.returncode:
        tail = "\n".join(log_path.read_text(encoding="utf-8", errors="replace").splitlines()[-30:])
        raise subprocess.CalledProcessError(completed.returncode, cmd, output=tail)


def _find_sra(run: str) -> Path:
    matches = list(SRA_DIR.glob(f"**/{run}.sra"))
    if matches:
        return matches[0]
    raise FileNotFoundError(f"prefetch terminou, mas não encontrei {run}.sra dentro de {SRA_DIR}")


def _gzip_fastq(path: Path, log_path: Path) -> Path:
    gz = Path(str(path) + ".gz")
    if gz.exists():
        return gz
    compressor = shutil.which("pigz") or shutil.which("gzip")
    if compressor is None:
        raise EnvironmentError("Instale gzip ou pigz para comprimir os FASTQs.")
    cmd = [compressor, "-f"]
    if Path(compressor).name.lower().startswith("pigz"):
        cmd += ["-p", str(THREADS)]
    cmd.append(str(path))
    _run(cmd, log_path)
    if not gz.is_file():
        raise FileNotFoundError(f"A compressão não produziu {gz}")
    return gz


def process_run(run: str) -> str:
    log_path = LOG_DIR / f"{run}.log"
    quant_file = QUANT_DIR / run / "quant.sf"
    if quant_file.is_file():
        logging.info("%s já quantificado; pulando.", run)
        return "already_quantified"

    run_raw = RAW_DIR / run
    run_trim = TRIM_DIR / run
    run_tmp = TMP_DIR / run
    run_raw.mkdir(parents=True, exist_ok=True)
    run_trim.mkdir(parents=True, exist_ok=True)
    run_tmp.mkdir(parents=True, exist_ok=True)
    SRA_DIR.mkdir(parents=True, exist_ok=True)
    QUANT_DIR.mkdir(parents=True, exist_ok=True)
    QC_DIR.mkdir(parents=True, exist_ok=True)

    try:
        sra = None
        try:
            sra = _find_sra(run)
        except FileNotFoundError:
            pass
        if sra is None:
            logging.info("Baixando %s", run)
            _run(["prefetch", run, "--output-directory", str(SRA_DIR)], log_path)
            sra = _find_sra(run)

        # Do not re-extract when FASTQ files from a prior interrupted run are present.
        fastqs = list(run_raw.glob("*.fastq")) + list(run_raw.glob("*.fastq.gz"))
        if not fastqs:
            _run(
                ["fasterq-dump", "--split-files", "--threads", str(THREADS), "--temp", str(run_tmp), "--outdir", str(run_raw), str(sra)],
                log_path,
            )
        # Keep fasterq-dump scratch files off the WSL /tmp tmpfs, which is often small.
        shutil.rmtree(run_tmp, ignore_errors=True)
        raw_r1 = run_raw / f"{run}_1.fastq"
        raw_r2 = run_raw / f"{run}_2.fastq"
        raw_single = run_raw / f"{run}.fastq"
        layout = "paired" if raw_r1.exists() or Path(str(raw_r1) + ".gz").exists() else "single"
        if layout == "paired":
            r1 = _gzip_fastq(raw_r1, log_path) if raw_r1.exists() else Path(str(raw_r1) + ".gz")
            r2 = _gzip_fastq(raw_r2, log_path) if raw_r2.exists() else Path(str(raw_r2) + ".gz")
            if not r1.is_file() or not r2.is_file():
                raise FileNotFoundError(f"FASTQ pareado incompleto para {run}: {r1}, {r2}")
            trim_cmd = ["trim_galore", "--paired", "--fastqc", "-j", str(THREADS), "--output_dir", str(run_trim), str(r1), str(r2)]
            trim_outputs = (run_trim / f"{run}_1_val_1.fq.gz", run_trim / f"{run}_2_val_2.fq.gz")
        else:
            if not raw_single.exists() and not Path(str(raw_single) + ".gz").exists():
                # fasterq-dump can emit one file with a spot suffix for unusual accessions.
                candidates = list(run_raw.glob("*.fastq")) + list(run_raw.glob("*.fastq.gz"))
                if len(candidates) != 1:
                    raise FileNotFoundError(f"Não consegui identificar FASTQ single-end para {run}: {candidates}")
                raw_single = candidates[0]
            single = _gzip_fastq(raw_single, log_path) if raw_single.exists() else Path(str(raw_single) + ".gz")
            trim_cmd = ["trim_galore", "--fastqc", "-j", str(THREADS), "--output_dir", str(run_trim), str(single)]
            trim_outputs = (run_trim / f"{single.name.removesuffix('.fastq.gz')}_trimmed.fq.gz",)
        _run(trim_cmd, log_path)
        if not all(p.is_file() for p in trim_outputs):
            raise FileNotFoundError(f"Trim Galore não gerou os FASTQs esperados: {trim_outputs}")

        salmon_cmd = ["salmon", "quant", "-i", str(REFERENCE_INDEX), "-l", "A", "-p", str(THREADS), "--validateMappings", "--seqBias", "--gcBias", "-o", str(QUANT_DIR / run)]
        if layout == "paired":
            salmon_cmd += ["-1", str(trim_outputs[0]), "-2", str(trim_outputs[1])]
        else:
            salmon_cmd += ["-r", str(trim_outputs[0])]
        _run(salmon_cmd, log_path)
        if not quant_file.is_file():
            raise FileNotFoundError(f"Salmon terminou sem criar {quant_file}")

        # Copy Trim Galore/FastQC reports into a stable per-run QC folder.
        reports = list(run_trim.glob("*report*.txt")) + list(run_trim.glob("*_fastqc.html")) + list(run_trim.glob("*_fastqc.zip"))
        run_qc = QC_DIR / run
        run_qc.mkdir(parents=True, exist_ok=True)
        for report in reports:
            shutil.copy2(report, run_qc / report.name)

        if not KEEP_INTERMEDIATES:
            shutil.rmtree(run_raw, ignore_errors=True)
            shutil.rmtree(run_trim, ignore_errors=True)
            shutil.rmtree(SRA_DIR / run, ignore_errors=True)
        return layout
    except Exception:
        logging.exception("Falha ao processar %s; arquivos e logs foram preservados para retomada.", run)
        raise


def run_pipeline() -> pd.DataFrame:
    _setup_logging()
    meta = _load_metadata()
    _require_tools()
    statuses = []
    for row in meta.itertuples(index=False):
        run = row.sample_id
        logging.info("Processando %s (%s / %s)", run, row.study_id, row.condition)
        try:
            layout = process_run(run)
            status = "already_quantified" if layout == "already_quantified" else "ok"
            if status == "already_quantified":
                layout = ""
            statuses.append({"sample_id": run, "condition": row.condition, "study_id": row.study_id, "library_layout": layout, "status": status})
        except Exception as exc:  # Continue with other runs; failures are summarized at the end.
            statuses.append({"sample_id": run, "condition": row.condition, "study_id": row.study_id, "library_layout": "", "status": f"failed: {exc}"})
    status_df = pd.DataFrame(statuses)
    status_df.to_csv(RUN_STATUS, index=False)
    enriched = meta.merge(status_df[["sample_id", "library_layout", "status"]], on="sample_id", how="left")
    enriched.to_csv(ENRICHED_METADATA, index=False)
    n_done = status_df.status.isin(["ok", "already_quantified"]).sum()
    logging.info("Pipeline concluído: %d/%d quantificados ou previamente prontos. Status: %s", n_done, len(status_df), RUN_STATUS)
    return enriched


if __name__ == "__main__":
    run_pipeline()
