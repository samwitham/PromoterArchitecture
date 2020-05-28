{
 "cells": [
  {
   "cell_type": "code",
   "execution_count": null,
   "metadata": {},
   "outputs": [],
   "source": [
    "#downloaded UCSC bigwig2wig from http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 1,
   "metadata": {},
   "outputs": [],
   "source": [
    "../../software/UCSC_utils/./bigWigToWig ../../data/ATAC-seq/potter2018/GSE116287_Roots_NaOH_Merged.bw ../../data/ATAC-seq/potter2018/Roots_NaOH_Merged.wig"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 2,
   "metadata": {},
   "outputs": [],
   "source": [
    "../../software/UCSC_utils/./bigWigToWig ../../data/ATAC-seq/potter2018/GSE116287_Shoots_NaOH_Merged.bw ../../data/ATAC-seq/potter2018/Shoots_NaOH_Merged.wig"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 3,
   "metadata": {},
   "outputs": [
    {
     "name": "stdout",
     "output_type": "stream",
     "text": [
      "\n",
      "Command 'convert2bed' not found, but can be installed with:\n",
      "\n",
      "sudo apt install bedops\n",
      "\n"
     ]
    },
    {
     "ename": "",
     "evalue": "127",
     "output_type": "error",
     "traceback": []
    }
   ],
   "source": [
    "convert2bed --do-not-sort -i wig < ../../data/ATAC-seq/potter2018/Roots_NaOH_Merged.wig > ../../data/ATAC-seq/potter2018/Roots_NaOH_Merged.bed"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 4,
   "metadata": {},
   "outputs": [
    {
     "name": "stdout",
     "output_type": "stream",
     "text": [
      "\n",
      "Command 'convert2bed' not found, but can be installed with:\n",
      "\n",
      "sudo apt install bedops\n",
      "\n"
     ]
    },
    {
     "ename": "",
     "evalue": "127",
     "output_type": "error",
     "traceback": []
    }
   ],
   "source": [
    "convert2bed --do-not-sort -i wig < ../../data/ATAC-seq/potter2018/GSE116287_Roots_NaOH_Merged.bw > ../../data/ATAC-seq/potter2018/Roots_NaOH_Merged.bed"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": null,
   "metadata": {},
   "outputs": [],
   "source": [
    "wig2bed < ../../data/ATAC-seq/potter2018/GSE116287_Roots_NaOH_Merged.bw > ../../data/ATAC-seq/potter2018/Roots_NaOH_Merged.bed"
   ]
  }
 ],
 "metadata": {
  "kernelspec": {
   "display_name": "Bash",
   "language": "bash",
   "name": "bash"
  },
  "language_info": {
   "codemirror_mode": "shell",
   "file_extension": ".sh",
   "mimetype": "text/x-sh",
   "name": "bash"
  }
 },
 "nbformat": 4,
 "nbformat_minor": 4
}
