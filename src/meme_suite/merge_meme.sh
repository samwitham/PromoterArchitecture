#!/bin/sh
## merges meme files together - specifically the DAP-seq meme files, using meme2meme from meme-suite
#$1 directory containing lots of sub directories containing .txt files containing motifs

meme2meme ${1}/**/*.txt > ../../data/FIMO/motif_data/${2}