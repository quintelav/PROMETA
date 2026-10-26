#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
PROMETA - Script de Expressão Diferencial e Metanálise (Método de Bono, 2021)
Focado estritamente no cálculo quantitativo das métricas de On-Score e On-Ratio.
Isolado de rotinas de anotação ou de leitura de ficheiros FASTA/GFF3.
"""

import os
import glob
import sys
import numpy as np
import pandas as pd

def carregar_dados_salmon(quant_dir):
    """
    Varre a pasta de resultados do Salmon e carrega os ficheiros quant.sf.
    Retorna um DataFrame unificado contendo as abundâncias TPM de todas as amostras.
    """
    if not os.path.exists(quant_dir):
        print(f"❌ [ERRO] A pasta '{quant_dir}' não existe no diretório atual!")
        print(f"Diretório atual de trabalho: {os.getcwd()}")
        print("Dica: Certifique-se de executar o script na raiz do projeto onde os dados do Salmon foram salvos.")
        return None

    dados_tpm = {}
    padrao_busca = os.path.join(quant_dir, "*", "quant.sf")
    arquivos = glob.glob(padrao_busca)
    
    if not arquivos:
        print(f"❌ [ERRO] Nenhum ficheiro quant.sf encontrado dentro de '{quant_dir}'!")
        print("Dica: Execute o pipeline de download/quantificação para as suas amostras primeiro.")
        return None
        
    print(f"📂 Encontrados {len(arquivos)} ficheiros de quantificação do Salmon.")
    
    for caminho in arquivos:
        sra_id = os.path.basename(os.path.dirname(caminho))
        try:
            df = pd.read_csv(caminho, sep='\t')
            if 'Name' in df.columns and 'TPM' in df.columns:
                # Usa o ID do transcrito como índice e extrai a coluna de TPM
                dados_tpm[sra_id] = df.set_index('Name')['TPM']
        except Exception as e:
            print(f"❌ Erro ao ler o ficheiro {caminho}: {e}")
            
    if not dados_tpm:
        print("❌ Não foi possível extrair dados válidos de TPM de nenhuma amostra.")
        return None
        
    df_expressao = pd.DataFrame(dados_tpm)
    return df_expressao

def calcular_on_ratio_on_score(df_tpm, metadados_csv, threshold_log2ratio=1.0):
    """
    Executa a Metanálise de Bono (2021) de forma estritamente matemática:
    1. Calcula Log2 Ratio por estudo: Log2((TPM_Tratado + 1) / (TPM_Control + 1))
    2. Calcula o On-Score (UP_Score / DOWN_Score) somando a ocorrência de atividade significativa.
    3. Calcula o Bono Net Score (UP_Score - DOWN_Score).
    4. Calcula o On-Ratio (UP_Score / Total de Estudos).
    """
    print(f"\n🧪 Iniciando análise quantitativa de On-Score e On-Ratio (Bono, 2021)...")
    print(f"🔹 Limiar de corte considerado: Log2(Ratio) >= {threshold_log2ratio} (ou <= -{threshold_log2ratio})")
    
    if not os.path.exists(metadados_csv):
        print(f"❌ [ERRO] O ficheiro de metadados '{metadados_csv}' não foi encontrado!")
        return None
        
    try:
        meta = pd.read_csv(metadados_csv)
    except Exception as e:
        print(f"❌ Erro ao ler o ficheiro de metadados {metadados_csv}: {e}")
        return None
        
    if 'study_id' not in meta.columns or 'sample_id' not in meta.columns or 'condition' not in meta.columns:
        print("❌ [ERRO] O ficheiro de metadados deve conter obrigatoriamente as colunas 'sample_id', 'condition' e 'study_id'!")
        return None
        
    estudos = meta['study_id'].unique()
    print(f"🧬 Estudos/Inseticidas detectados nos metadados: {list(estudos)}")
    
    total_estudos = len(estudos)
    colunas_ratios = []
    
    # Inicializa DataFrame que guardará os resultados da metanálise
    df_bono = pd.DataFrame(index=df_tpm.index)
    
    for estudo in estudos:
        meta_sub = meta[meta['study_id'] == estudo]
        
        controles = meta_sub[meta_sub['condition'].str.lower() == 'control']['sample_id'].tolist()
        tratados = meta_sub[meta_sub['condition'].str.lower() == 'treatment']['sample_id'].tolist()
        
        # Filtrar apenas as amostras do metadados que realmente existem nos resultados do Salmon
        controles = [c for c in controles if c in df_tpm.columns]
        tratados = [t for t in tratados if t in df_tpm.columns]
        
        if not controles or not tratados:
            print(f"⚠️  [AVISO] Amostras insuficientes mapeadas no Salmon para o estudo '{estudo}'. Pulando.")
            continue
            
        print(f"   📊 Processando '{estudo}': {len(controles)} controles vs {len(tratados)} tratados.")
        
        # Médias de TPM de cada grupo
        tpm_control = df_tpm[controles].mean(axis=1)
        tpm_treatment = df_tpm[tratados].mean(axis=1)
        
        # Fórmula matemática de Bono (2021) com pseudocount de +1.0 para suavização de baixas contagens
        log2_ratio = np.log2((tpm_treatment + 1.0) / (tpm_control + 1.0))
        
        # Salva as abundâncias médias e o log2 ratio correspondente na tabela final
        df_bono[f'TPM_Control_{estudo}'] = tpm_control
        df_bono[f'TPM_Treatment_{estudo}'] = tpm_treatment
        df_bono[f'Log2Ratio_{estudo}'] = log2_ratio
        colunas_ratios.append(f'Log2Ratio_{estudo}')

    if not colunas_ratios:
        print("❌ [ERRO] Falha crítica: Nenhum estudo pôde ser processado. Verifique os IDs de amostras nos metadados.")
        return None

    # --- CÁLCULO DAS MÉTRICAS DE CONSISTÊNCIA DE BONO ---
    print("\n📈 Computando On-Scores e On-Ratios ao longo das comparações...")
    
    # Inicializa as matrizes auxiliares de somatório
    up_score_matrix = np.zeros(len(df_bono))
    down_score_matrix = np.zeros(len(df_bono))
    
    for col in colunas_ratios:
        # Atribui +1 para cada comparação em que o transcrito ultrapassa o limiar estipulado (ativo)
        up_score_matrix += (df_bono[col] >= threshold_log2ratio).astype(int)
        down_score_matrix += (df_bono[col] <= -threshold_log2ratio).astype(int)
        
    df_bono['Bono_UP_Score'] = up_score_matrix.astype(int)
    df_bono['Bono_DOWN_Score'] = down_score_matrix.astype(int)
    
    # Bono Net Score: Indica a tendência resultante final do gene (Positivo = Ativado / Negativo = Reprimido)
    df_bono['Bono_Net_Score'] = df_bono['Bono_UP_Score'] - df_bono['Bono_DOWN_Score']
    
    # On-Ratio: Representa a frequência relativa de atividade do gene entre todos os ensaios analisados
    df_bono['On_Ratio_UP'] = df_bono['Bono_UP_Score'] / total_estudos
    df_bono['On_Ratio_DOWN'] = df_bono['Bono_DOWN_Score'] / total_estudos
    
    return df_bono

if __name__ == "__main__":
    print("==================================================")
    print("   PROMETA - Pipeline Metanálise Bono (2021)       ")
    print("   ETAPA 1: EXPRESSÃO DIFERENCIAL (GFF3-FREE)     ")
    print("==================================================")
    
    diretorio_quant = "quant"
    metadados_geral = "metadata_geral.csv"
    
    # Executa o carregamento das abundâncias calculadas pelo Salmon
    matriz_tpm = carregar_dados_salmon(diretorio_quant)
    
    if matriz_tpm is not None:
        # Garante a criação da pasta quant se ela ainda não existir
        os.makedirs(diretorio_quant, exist_ok=True)
        
        # Salva a matriz bruta de TPM consolidada para backup e auditoria do projeto
        matriz_tpm.to_csv(f"{diretorio_quant}/matriz_tpm_geral.csv")
        print(f"💾 Matriz de TPM geral unificada guardada em '{diretorio_quant}/matriz_tpm_geral.csv'.")
        
        # Calcula as métricas matemáticas On-Score e On-Ratio da metanálise
        # threshold_log2ratio=1.0 representa uma variação mínima de 2 vezes (Fold-Change >= 2)
        df_metanalise_quant = calcular_on_ratio_on_score(
            matriz_tpm, metadados_geral, threshold_log2ratio=1.0
        )
        
        if df_metanalise_quant is not None:
            # Salva os resultados quantitativos puros da metanálise
            output_metanalise = f"{diretorio_quant}/metanalise_completa_bono.csv"
            df_metanalise_quant.to_csv(output_metanalise)
            print(f"💾 Resultados consolidados salvos em '{output_metanalise}'.")
            
            # --- FILTRAGEM DOS CANDIDATOS ATIVOS GERAIS ---
            # Considera ativo qualquer transcrito que alterou sua expressão em pelo menos 1 comparação
            score_corte = 1
            
            up_geral = df_metanalise_quant[df_metanalise_quant['Bono_UP_Score'] >= score_corte]
            down_geral = df_metanalise_quant[df_metanalise_quant['Bono_DOWN_Score'] >= score_corte]
            
            # Ordena os transcritos pelo On-Score de forma decrescente (destacando os mais consistentes no topo)
            up_geral = up_geral.sort_values(by='Bono_UP_Score', ascending=False)
            down_geral = down_geral.sort_values(by='Bono_DOWN_Score', ascending=False)
            
            up_geral.to_csv(f"{diretorio_quant}/candidatos_up_onscore_geral.csv")
            down_geral.to_csv(f"{diretorio_quant}/candidatos_down_onscore_geral.csv")
            
            print("\n==================================================")
            print("📊 ESTATÍSTICAS QUANTITATIVAS DA METANÁLISE:")
            print("==================================================")
            print(f"   📈 Transcritos Ativos Totais (Super-expressos): {len(up_geral)}")
            print(f"   📉 Transcritos Ativos Totais (Sub-expressos): {len(down_geral)}")
            
            # Detalhamento de frequência de consistência molecular
            total_estudos_detetados = len(df_metanalise_quant.columns) // 3
            for s in range(1, total_estudos_detetados + 1):
                total_s_up = len(df_metanalise_quant[df_metanalise_quant['Bono_UP_Score'] == s])
                total_s_down = len(df_metanalise_quant[df_metanalise_quant['Bono_DOWN_Score'] == s])
                if total_s_up > 0 or total_s_down > 0:
                    print(f"   🔹 Ativo em exatamente {s} estudo(s) -> Super-expressos: {total_s_up} | Sub-expressos: {total_s_down}")
            
            print("\n🎉 Etapa 1 Concluída! Suas matrizes quantitativas de expressão diferencial estão prontas.")
            print("Próximo passo: Use o script de anotação funcional 'anotar_resultados.py' para classificar estes transcritos!")
            
        else:
            print("❌ Erro durante o processamento das comparações estatísticas.")
    else:
        print("❌ Falha crítica ao ler abundâncias. Certifique-se de que os dados do Salmon estão na pasta 'quant/'.")
