#!/bin/bash
set -e  # Interrompe o script se qualquer comando falhar

# 1. Variáveis de Ambiente e Diretórios Absolutos
PROJECT_ROOT=$(pwd)
KAGGLE_USERNAME="calistu"

NOTEBOOK_ID_BASE="notebook-annexyn-a11"
NOTEBOOK_ID="calistu/$NOTEBOOK_ID_BASE"
TIMESTAMP=$(date +%s)
DATASET_SLUG="anexina-frames-$TIMESTAMP"

# Usando o PROJECT_ROOT para garantir caminhos absolutos precisos
EXPORT_DIR="$PROJECT_ROOT/results/$(date +'%Y-%m-%d')/exported_frames"
KAGGLE_DIR="$PROJECT_ROOT/kaggle"

# ATENÇÃO: Coloque o nome real do seu arquivo .ipynb aqui

NOTEBOOK_FILE="$NOTEBOOK_ID_BASE.ipynb" 

# 2. Exportação dos Frames (Bioemu / Anexina A11)
echo "--- Iniciando export_frames.py ---"
rm -rf $EXPORT_DIR
python3 scripts/export_frames.py

# 3. Preparação e Criação do Dataset (Com Compressão)
echo "--- Preparando Dataset: $DATASET_SLUG ---"

if [ ! -d "$EXPORT_DIR" ]; then
  echo "Erro: O diretório $EXPORT_DIR não foi encontrado."
  exit 1
fi

# Cria uma pasta temporária apenas para o upload
STAGING_DIR="$PROJECT_ROOT/dataset_staging_$TIMESTAMP"
mkdir -p "$STAGING_DIR"

# Compacta os frames gerados pelo bioemu direto para a pasta de staging
echo "Compactando arquivos de $EXPORT_DIR..."
cd "$EXPORT_DIR"
zip -r -q "$STAGING_DIR/bioemu_frames.zip" .

# Vai para a pasta de staging gerar o metadata e fazer o upload
cd "$STAGING_DIR"

cat <<EOF > dataset-metadata.json
{
  "title": "$DATASET_SLUG",
  "id": "$KAGGLE_USERNAME/$DATASET_SLUG",
  "licenses": [{"name": "CC0-1.0"}]
}
EOF

# Faz o upload apenas do .zip e do json
echo "Enviando para o Kaggle..."
kaggle datasets create -p . -u

# Limpa o diretório temporário para não acumular lixo
cd "$PROJECT_ROOT"
rm -rf "$STAGING_DIR"

# 4. Configuração e Execução do Kernel
echo "--- Configurando Kernel: $NOTEBOOK_ID ---"

mkdir -p $KAGGLE_DIR
cd "$KAGGLE_DIR"

cat <<EOF > kernel-metadata.json
{
  "id": "$NOTEBOOK_ID",
  "title": "Notebook Annexyn A11",
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

kaggle kernels pull calistu/notebook-anexina-a11

# Aumentado para 30s para garantir que o Kaggle indexe o dataset
echo "Aguardando 30 segundos para processamento do dataset no backend do Kaggle..."
sleep 30

echo "--- Disparando execução no Kaggle ---"
kaggle kernels push -p .

echo "--- Sucesso. Acompanhe em: https://www.kaggle.com/$NOTEBOOK_ID ---"