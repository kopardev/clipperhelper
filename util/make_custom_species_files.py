import gtfparse
import pandas
import sys
import argparse

parser = argparse.ArgumentParser(description='Create clipper input files for custom species')
parser.add_argument('--gtf', dest='gtf', type=str, required=True,
                    help='input GTF file')
parser.add_argument('--species', dest='species', type=str, required=True,
                    help='species name to be used as output file prefix')
args = parser.parse_args()
debug=0

gene_outfile="%s_genes.bed"%(args.species)
exon_outfile="%s_exons.bed"%(args.species)
gtf_outfile="%s.AS.STRUCTURE.COMPILED.gff"%(args.species)

df = gtfparse.read_gtf(args.gtf)
seqnames = df.seqname.unique()
bigdict = dict()
if debug==1: print("list of genes=",list(df.gene_id.unique()))

gene2loc=dict()
genedata=[]
exondata=[]

for seq in seqnames:
    print(seq)
    df_seq=df[df['seqname']==seq]
    for gene in list(df_seq.gene_id.unique()):
        y=df_seq[(df_seq["gene_id"]==gene) & (df_seq["feature"]=="gene")]
        gene_chrom=y['seqname'].iloc[0]
        gene_start=str(y['start'].iloc[0])
        gene_end=str(y['end'].iloc[0])
        gene_strand=y['strand'].iloc[0]
        gene2loc[gene]="##".join([gene_chrom,gene_start,gene_end,gene_strand])
        gene_start=str(int(gene_start)-1)
        genedata.append([gene_chrom,gene_start,gene_end,gene,"0",gene_strand])
        bigdict[gene]=list()
        x=df_seq[(df_seq["gene_id"]==gene) & (df_seq["feature"]!="gene")]
        if debug==1: print("gene=",gene)
        if debug==1: print("list of transcripts=",list(x.transcript_id.unique()))
        for transcript in list(x.transcript_id.unique()):
            y=x[(x["transcript_id"]==transcript) & (x["feature"]=="exon")]
            l=0
            for i,exonrow in y.iterrows():
                exon_start=str(exonrow["start"]-1)
                exon_end=str(exonrow["end"])
                exondata.append([gene_chrom,exon_start,exon_end,gene,"0",gene_strand])
                l+=(exonrow["end"]-(exonrow["start"]-1))
            if debug==1: print("transcript=",transcript)
            if debug==1: print("l=",l)
            bigdict[gene].append(l)
        if debug==1: print(bigdict)

genedf=pandas.DataFrame(genedata,columns=['chrom','start','end','geneid','score','strand'])
genedf=genedf.astype({'start':int,'end':int})
genedf=genedf.sort_values(by=['chrom','start','end'])
genedf.to_csv(gene_outfile,sep="\t",header=False,index=False)

exondf=pandas.DataFrame(exondata,columns=['chrom','start','end','geneid','score','strand'])
exondf=exondf.drop_duplicates()
exondf=exondf.astype({'start':int,'end':int})
exondf=exondf.sort_values(by=['chrom','start','end'])
exondf.to_csv(exon_outfile,sep="\t",header=False,index=False)

gtfdata=[]
for gene in list(df.gene_id.unique()):
    maxlen=max(bigdict[gene])
    chrom,start,end,strand=gene2loc[gene].split("##")
    premrna_length=int(end)-(int(start)-1)
    attribute_text="gene_id=%s;mrna_length=%d;premrna_length=%d"%(gene,maxlen,premrna_length)
    r=[chrom,"AS_STRUCTURE","gene",start,end,".",strand,".",attribute_text]
    gtfdata.append(r)

gtfdf=pandas.DataFrame(gtfdata)
gtfdf.to_csv(gtf_outfile,sep="\t",header=False,index=False)