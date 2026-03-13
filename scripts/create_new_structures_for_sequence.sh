#!/bin/bash
set -e  # Interrompe o script se qualquer comando falhar

# 1. Variáveis de Ambiente e Diretórios Absolutos
PROJECT_ROOT=$(pwd)
KAGGLE_USERNAME="calistu"

NOTEBOOK_ID_BASE="notebook-annexyn-a11-bioemu"
NOTEBOOK_ID="calistu/$NOTEBOOK_ID_BASE"
TIMESTAMP=$(date +%s)
DATASET_SLUG="anexina-frames-$TIMESTAMP"

# Usando o PROJECT_ROOT para garantir caminhos absolutos precisos
EXPORT_DIR="$PROJECT_ROOT/results/$(date +'%Y-%m-%d')/bioemu_samples"
KAGGLE_DIR="$PROJECT_ROOT/kaggle"

# ATENÇÃO: Coloque o nome real do seu arquivo .ipynb aqui

NOTEBOOK_FILE="$NOTEBOOK_ID_BASE.ipynb" 

echo "--- Configurando Kernel: $NOTEBOOK_ID ---"

mkdir -p $KAGGLE_DIR
cd "$KAGGLE_DIR"

cat <<EOF > kernel-metadata.json
{
  "id": "$NOTEBOOK_ID",
  "title": "Notebook Annexyn A11 Bioemu",
  "code_file": "$NOTEBOOK_FILE",
  "language": "python",
  "kernel_type": "notebook",
  "is_private": "true",
  "enable_gpu": "true",
  "enable_tpu": "false",
  "enable_internet": "true",
  "dataset_sources": ["$KAGGLE_USERNAME/$DATASET_SLUG"],
  "kernel_sources": [],
  "competition_sources": []
}
EOF

# kaggle kernels pull calistu/notebook-anexina-a11-bioemu
#kaggle kernels init

echo "--- Disparando execução no Kaggle ---"
kaggle kernels push -p .

echo "--- Sucesso. Acompanhe em: https://www.kaggle.com/$NOTEBOOK_ID ---"


echo "--- Aguardando a execução do kernel no Kaggle ---"

while true; do
  # Captura o status atual do kernel
  STATUS=$(kaggle kernels status $NOTEBOOK_ID)
  
  # Verifica se a palavra "complete" está no retorno
  if echo "$STATUS" | grep -qi "complete"; then
    echo "Execução finalizada com sucesso."
    break
  fi
  
  # Verifica se deu erro para não ficar em loop infinito
  if echo "$STATUS" | grep -qi "error\|cancel"; then
    echo "A execução falhou ou foi cancelada. Verifique no Kaggle."
    exit 1
  fi
  
  echo "Kernel ainda rodando. Checando novamente em 60 segundos..."
  sleep 60
done


mkdir -p $EXPORT_DIR 
kaggle kernels output $KAGGLE_USERNAME/$DATASET_SLUG -p $EXPORT_DIR
    