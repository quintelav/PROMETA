#!/usr/bin/env bash
set -Eeuo pipefail
export PATH="$HOME/.local/bin:$PATH"

# Run from WSL. The analysis and large sequencing files live in Linux's filesystem,
# while this script copies the current project sources from the Windows workspace.
WINDOWS_PROJECT="${PROMETA_WINDOWS_PROJECT:-/mnt/c/Users/quint/OneDrive/Documentos/ChatGPT/PROMETA}"
PROJECT_DIR="${PROMETA_PROJECT_DIR:-$HOME/projects/PROMETA}"
LOCAL_METADATA="$WINDOWS_PROJECT/metadata_geral.csv"
WSL_METADATA="$HOME/projects/tuta_rnaseq/metadata/metadata_geral.csv"

if [[ ! -d "$WINDOWS_PROJECT" ]]; then
  echo "Pasta do projeto Windows não encontrada: $WINDOWS_PROJECT" >&2
  echo "Defina PROMETA_WINDOWS_PROJECT com o caminho correto em /mnt/<unidade>/..." >&2
  exit 1
fi

if command -v micromamba >/dev/null 2>&1; then
  ENV_MANAGER="micromamba"
elif command -v mamba >/dev/null 2>&1; then
  ENV_MANAGER="mamba"
elif command -v conda >/dev/null 2>&1; then
  ENV_MANAGER="conda"
else
  echo "Não encontrei conda, mamba ou micromamba no WSL." >&2
  echo "Instale Miniforge/Mambaforge ou micromamba e rode este script novamente." >&2
  exit 1
fi

mkdir -p "$PROJECT_DIR" "$PROJECT_DIR/reference"
for file in projeto_paralelo_vitor_quintela.ipynb prometa_pipeline.py analise_diferencial.R analise_meta_bono.py environment.yml setup_wsl.sh README.md README_WSL.md; do
  if [[ -f "$WINDOWS_PROJECT/$file" ]]; then
    cp -f "$WINDOWS_PROJECT/$file" "$PROJECT_DIR/$file"
  else
    echo "Arquivo do projeto ausente: $WINDOWS_PROJECT/$file" >&2
    exit 1
  fi
done
if [[ -f "$WINDOWS_PROJECT/.vscode/extensions.json" ]]; then
  mkdir -p "$PROJECT_DIR/.vscode"
  cp -f "$WINDOWS_PROJECT/.vscode/extensions.json" "$PROJECT_DIR/.vscode/extensions.json"
fi

if [[ -f "$LOCAL_METADATA" ]]; then
  cp -f "$LOCAL_METADATA" "$PROJECT_DIR/metadata_geral.csv"
elif [[ -f "$WSL_METADATA" ]]; then
  cp -f "$WSL_METADATA" "$PROJECT_DIR/metadata_geral.csv"
elif [[ ! -f "$PROJECT_DIR/metadata_geral.csv" ]]; then
  echo "Não encontrei metadata_geral.csv em $LOCAL_METADATA nem em $WSL_METADATA." >&2
  exit 1
fi

# Reuse the existing Tuta reference without copying the large Salmon index.
REFERENCE_SOURCE="$HOME/projects/tuta_rnaseq/reference"
if [[ ! -e "$PROJECT_DIR/reference/salmon_index" && -d "$REFERENCE_SOURCE/salmon_index" ]]; then
  ln -s "$REFERENCE_SOURCE/salmon_index" "$PROJECT_DIR/reference/salmon_index"
fi
if [[ ! -e "$PROJECT_DIR/reference/tuta_transcripts.fna" && -f "$REFERENCE_SOURCE/tuta_transcripts.fna" ]]; then
  ln -s "$REFERENCE_SOURCE/tuta_transcripts.fna" "$PROJECT_DIR/reference/tuta_transcripts.fna"
fi

if "$ENV_MANAGER" env list | awk '{print $1}' | grep -qx 'prometa-rnaseq'; then
  "$ENV_MANAGER" env update -n prometa-rnaseq -f "$PROJECT_DIR/environment.yml"
else
  "$ENV_MANAGER" env create -f "$PROJECT_DIR/environment.yml"
fi

"$ENV_MANAGER" run -n prometa-rnaseq python -m ipykernel install --user \
  --name prometa-rnaseq --display-name "Python (PROMETA RNA-seq)"

if [[ -f "$PROJECT_DIR/reference/tuta_transcripts.fna" && ! -f "$PROJECT_DIR/reference/tx2gene.csv" ]]; then
  "$ENV_MANAGER" run -n prometa-rnaseq python - "$PROJECT_DIR/reference/tuta_transcripts.fna" \
    "$PROJECT_DIR/reference/tx2gene.csv" <<'PY'
import csv
import re
import sys

fasta_path, output_path = sys.argv[1:]
pattern = re.compile(r"rna-gnl\|WGS:[^|]+\|(Tabs_g[0-9]+)\.t[0-9]+$")
rows = []
with open(fasta_path, encoding="utf-8") as fasta:
    for line in fasta:
        if line.startswith(">"):
            transcript = line[1:].split()[0]
            match = pattern.fullmatch(transcript)
            if match is None:
                raise SystemExit(f"Cabeçalho fora do padrão esperado no FASTA: {transcript}")
            rows.append((transcript, match.group(1)))
if not rows:
    raise SystemExit("Nenhum transcrito encontrado no FASTA.")
with open(output_path, "w", newline="", encoding="utf-8") as output:
    writer = csv.writer(output)
    writer.writerow(["TXNAME", "GENEID"])
    writer.writerows(rows)
print(f"Mapa TXNAME/GENEID criado com {len(rows)} transcritos: {output_path}")
PY
fi

mkdir -p "$PROJECT_DIR/logs" "$PROJECT_DIR/qc" "$PROJECT_DIR/quant" \
  "$PROJECT_DIR/raw_data" "$PROJECT_DIR/sra_bruto"

cat <<EOF

Ambiente PROMETA preparado.
Projeto: $PROJECT_DIR
Gerenciador: $ENV_MANAGER

Próximos passos:
  1. Confirme se o índice Salmon e tx2gene.csv estão em $PROJECT_DIR/reference.
     O setup reutiliza o índice Tuta de ~/projects/tuta_rnaseq/reference quando disponível.
  2. Abra no VS Code com: cd "$PROJECT_DIR" && code .
  3. No notebook, selecione o kernel "Python (PROMETA RNA-seq)".
  4. Rode primeiro a célula de download/QC/Salmon.

Para ativar o ambiente no terminal WSL: $ENV_MANAGER activate prometa-rnaseq
EOF

