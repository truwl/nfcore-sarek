set -e
echo "trying sarek"
export BD_ENV='stagingadmin'
python bdcli.py workflows generate_nf2wdl ./src/lib/nf2wdl/workflows/nfcore-sarek-5d0a0b97551a9700b9532ed4d32c4e9c543a998e/nfcore_sarek_schema.json sarek 5d0a0b97551a9700b9532ed4d32c4e9c543a998e
mkdir -p ../capanno/workflows/nfcore/sarek/5d0a0b97551a9700b9532ed4d32c4e9c543a998e/
cp ./src/lib/nf2wdl/workflows/nfcore-sarek-5d0a0b97551a9700b9532ed4d32c4e9c543a998e/inputs.json ./src/lib/nf2wdl/workflows/nfcore-sarek-5d0a0b97551a9700b9532ed4d32c4e9c543a998e/inputs-metadata.yaml ./src/lib/nf2wdl/workflows/nfcore-sarek-5d0a0b97551a9700b9532ed4d32c4e9c543a998e/sarek.wdl ~/Documents/dev/truwl/capanno/workflows/nfcore/sarek/5d0a0b97551a9700b9532ed4d32c4e9c543a998e/
python bdcli.py workflows add nfcore sarek 5d0a0b97551a9700b9532ed4d32c4e9c543a998e
git -C ./src/lib/nf2wdl/workflows/nfcore-sarek-5d0a0b97551a9700b9532ed4d32c4e9c543a998e/ add inputs.json *metadata.yaml sarek.wdl
git -C ./src/lib/nf2wdl/workflows/nfcore-sarek-5d0a0b97551a9700b9532ed4d32c4e9c543a998e/ commit -m 'inputs.json inputs-metadata.yaml sarek.wdl' -a
git -C ./src/lib/nf2wdl/workflows/nfcore-sarek-5d0a0b97551a9700b9532ed4d32c4e9c543a998e/ push origin master
echo "sarek success"
        