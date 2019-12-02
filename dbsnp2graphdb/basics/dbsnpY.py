import urllib.request
import os
import json
import bz2
import re
from dataclasses import dataclass
from typing import List
import itertools

""" 
@author: Lauren DeLong

The following script downloads dbSNP data for Chromosome Y. 
The script has been adapted from the bio2bel package and is intended to be used for a graph database tutorial 
with OrientDB. If you use this code please give credit to the bio2bel creators. 
If you want the full versions which parse all chromosomes and/or write data to a relational database 
please refer to: https://github.com/bio2bel/dbsnp 
"""


@dataclass
class SnpInfo:
    dbsnp_id: str
    assembly_id: str
    gen_type: str
    gene_name: str
    gene_abbr: str
    gene_id: str
    dna_change: str
    rnas: str
    rna_type: str
    rna_change: List[str]
    proteins: str
    prot_type: str
    aa_change: List[str]


def createGeneTable(cl, out, count):
    print(count,
          cl.gene_id,
          cl.gene_abbr,
          cl.gene_name,
          sep=',',
          file=out
          )


def createProteinTable(cl, out, count):
    print(count,
          cl.proteins,
          cl.prot_type,
          cl.gene_id,
          sep=',',
          file=out
          )


def createRNATable(cl, out, count):
    print(count,
          cl.rnas,
          cl.rna_type,
          cl.gene_id,
          cl.proteins,
          sep=',',
          file=out
          )


def main():
    # Here we begin the downloading of JSON files from the dbSNP database:

    url = 'ftp://ftp.ncbi.nih.gov/snp/latest_release/JSON/refsnp-chrY.json.bz2'
    path = # file name here
    if not os.path.exists(path):
        print('Beginning file download with urllib2...')
        urllib.request.urlretrieve(url, path)
        print('...Finished file download with urllib2.')

    # Here we parse through the files:
    print('Now decompressing and reading JSON.bz2 files with *bz2* and *json* ...')
    with bz2.BZ2File(path, 'rb') as f_in, \
            open('FileName/Gene-chrY.csv', 'w') as Gene, \
            open('/FileName/RNA-chrY.csv', 'w') as RNA, \
            open('FileName/Protein-chrY.csv', 'w') as Protein:

        print('id', 'entrez_id', 'symbol', 'name', sep=',', file=Gene)
        print('id', 'accession', 'type', 'gene_id', 'protein_id', sep=',', file=RNA)
        print('id', 'accession', 'type', 'gene_id', sep=',', file=Protein)

        gene_tuple_list = []
        protein_tuple_list = []
        rna_tuple_list = []

        for line in f_in:
            rs_obj = json.loads(line.decode('utf-8'))
            dbsnp_id = rs_obj['refsnp_id']  # the dbsnp id

            all_ann_list_raw = rs_obj['primary_snapshot_data'][
                'allele_annotations']  # these are the assembly annotations

            if len(all_ann_list_raw) >= 2:  # if it has sufficient info
                assembl_ann_list_raw = all_ann_list_raw[1]['assembly_annotation']  # against each assembly
                if len(assembl_ann_list_raw) != 0:  # if it contains gene info
                    gene_list_raw = assembl_ann_list_raw[0][
                        'genes']  # and each of the genes affected within each assembly
                    if len(gene_list_raw) > 0:
                        # Here I start extracting gene info:
                        for x, y, z in itertools.product(range(len(all_ann_list_raw)),
                                                         range(len(assembl_ann_list_raw)),
                                                         range(len(gene_list_raw))):
                            assembly_id = all_ann_list_raw[x]['assembly_annotation'][y]['seq_id']
                            if assembly_id[0:2] == 'AC':
                                gen_type = 'Complete genomic molecule, usually alternate assembly'
                            elif assembly_id[0:2] == 'NC':
                                gen_type = 'Complete genomic molecule, usually reference assembly'
                            elif assembly_id[0:2] == 'NG':
                                gen_type = 'Incomplete genomic region'
                            elif assembly_id[0:2] == 'NT':
                                gen_type = 'Contig or scaffold, clone-based or WGSa'
                            elif assembly_id[0:2] == 'NW':
                                gen_type = 'Contig or scaffold, primarily WGSa'
                            elif assembly_id[0:2] == 'NZ':
                                gen_type = 'Complete genomes and unfinished WGS data'
                            else:
                                gen_type = ''

                            gene_name = all_ann_list_raw[x]['assembly_annotation'][y]['genes'][z]['name']
                            gene_abbr = all_ann_list_raw[x]['assembly_annotation'][y]['genes'][z]['locus']
                            gene_id = all_ann_list_raw[x]['assembly_annotation'][y]['genes'][z]['id']
                            rna_list_raw = all_ann_list_raw[x]['assembly_annotation'][y]['genes'][z]['rnas']

                            gene_tuple = (gene_name, gene_abbr, gene_id)
                            if gene_tuple not in gene_tuple_list:
                                gene_tuple_list.append(gene_tuple)

                            for nuc in rna_list_raw:
                                if 'id' in nuc:
                                    rnas = nuc['id']  # the rna transcript affected by the mutation
                                    if rnas[0:2] == 'NM':
                                        rna_type = 'Protein-coding transcripts (usually curated)'
                                    elif rnas[0:2] == 'NR':
                                        rna_type = 'Non-protein-coding transcripts'
                                    elif rnas[0:2] == 'XM':
                                        rna_type = 'Predicted model protein-coding transcript'
                                    elif rnas[0:2] == 'XR':
                                        rna_type = 'Predicted model non-protein-coding transcript'
                                    else:
                                        rna_type = ''
                                else:
                                    rnas = ''
                                    rna_type = ''
                                if 'product_id' in nuc:
                                    proteins = nuc['product_id']  # the protein affected by the mutation
                                    if proteins[0:2] == 'AP':
                                        prot_type = 'Annotated on AC_ alternate assembly'
                                    elif proteins[0:2] == 'NP':
                                        prot_type = 'Associated with an NM_ or NC_ accession'
                                    elif proteins[0:2] == 'YP':
                                        prot_type = 'Annotated on genomic molecules without an instantiated ' \
                                                    'transcript record '
                                    elif proteins[0:2] == 'XP':
                                        prot_type = 'Predicted model, associated with an XM_ accession'
                                    elif proteins[0:2] == 'WP':
                                        prot_type = 'Non-redundant across multiple strains and species'
                                    else:
                                        prot_type = ''
                                else:
                                    proteins = ''
                                    prot_type = ''

                                protein_tuple = (gene_id, proteins, prot_type)
                                if (proteins != '') and (protein_tuple not in protein_tuple_list):
                                    protein_tuple_list.append(protein_tuple)
                                rna_tuple = (gene_id, rnas, rna_type, proteins)
                                if (rnas != '') and (rna_tuple not in rna_tuple_list):
                                    rna_tuple_list.append(rna_tuple)

                                # Here I parse through each hgvs entry and assign it to either a nuc. change or a.a. change
                                hgvs_entries = rs_obj['primary_snapshot_data']['placements_with_allele']
                                dna_change = []
                                aa_change = []
                                rna_change = []
                                for entry in hgvs_entries:
                                    for variant in entry['alleles']:
                                        hgvs = re.split(":[cgmnopr].", variant['hgvs'])
                                        if len(hgvs) > 1:
                                            if hgvs[0] == assembly_id:
                                                dna_change.append(hgvs[1])
                                            elif hgvs[0] == proteins:
                                                aa_change.append(hgvs[1])
                                            elif hgvs[0] == rnas:
                                                rna_change.append(hgvs[1])
                                            else:
                                                continue

                                max_list = max(len(dna_change), len(rna_change), len(aa_change))
                                for i in [dna_change, rna_change, aa_change]:
                                    diff = abs(max_list - len(i))
                                    i.extend(list(itertools.repeat('', diff)))

        for count, gene in enumerate(gene_tuple_list):
            snp_infos = SnpInfo('',
                                '',
                                '',
                                gene[0],
                                gene[1],
                                gene[2],
                                '',
                                '',
                                '',
                                [''],
                                '',
                                '',
                                ['']
                                )

            createGeneTable(snp_infos, Gene, count)

        for count, rna in enumerate(rna_tuple_list):
            snp_infos = SnpInfo('',
                                '',
                                '',
                                '',
                                '',
                                rna[0],
                                '',
                                rna[1],
                                rna[2],
                                [''],
                                rna[3],
                                '',
                                ['']
                                )

            createRNATable(snp_infos, RNA, count)

        for count, prot in enumerate(protein_tuple_list):
            snp_infos = SnpInfo('',
                                '',
                                '',
                                '',
                                '',
                                prot[0],
                                '',
                                '',
                                '',
                                [''],
                                prot[1],
                                prot[2],
                                ['']
                                )

            createProteinTable(snp_infos, Protein, count)

    print("Finished writing files to CSV.")


if __name__ == '__main__':
    main()
