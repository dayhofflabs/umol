#!/bin/sh
# Fetch the VEHICLe ring-system corpus into the local staging directory.
#
# Provenance: Pitt et al., "Heteroaromatic Rings of the Future",
# J. Med. Chem. 2009, 52, 9, 2952-2963, doi:10.1021/jm801513z.
# Source of record: ChEMBL FTP, pub/databases/chembl/VEHICLe/ (dated
# 2010-04-12). The FTP README carries no licence statement; the file is a
# local development instrument only and must never be committed or
# redistributed (discussion/196, "Corpus"). The staging directory lies under
# materials/, which is covered by the repository .gitignore.
set -eu

staging="$(dirname "$0")/../materials/aromaticity/vehicle"
base="https://ftp.ebi.ac.uk/pub/databases/chembl/VEHICLe"

mkdir -p "$staging"
curl -fsS "$base/README" -o "$staging/README"
curl -fsS "$base/VEHICLe.csv" -o "$staging/VEHICLe.csv"
wc -l "$staging/VEHICLe.csv"
