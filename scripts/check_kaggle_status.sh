NOTEBOOK_ID="notebook-annexyn-a11-bioemu"

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
