#!/usr/bin/env Rscript

args = commandArgs(trailingOnly=TRUE)

mosdepth <- args[1]
ref_gtf <- args[2]
stringtie_gtf <- args[3]
output <- args[4]

library(dplyr)
library(readr)
library(ggplot2)
library(cowplot)
library(tidyr)
library(ggthemes)


cov = read_tsv(mosdepth,
               col_names = c("chr", "start", "end", "coverage"))
genes = read_tsv(stringtie_gtf,
                 col_names = c("chr", "start", "end", "strand", "exons",
                               "tx", "gene", "name", "tpm")) %>%
    filter(tpm != ".") %>% 
    mutate(name=ifelse(name==".",paste0("tx", 1:n()),name)) %>% 
    mutate(ntx = 1:n(),
           direction=ifelse(strand=="+", "last", "first")) 
ref_genes = read_tsv(ref_gtf,
                 col_names = c("chr", "start", "end", "strand", "exons",
                                "name")) %>% 
    mutate(name=ifelse(name==".",paste0("tx", 1:n()),name)) %>% 
    mutate(ntx = 1:n(),
           direction=ifelse(strand=="+", "last", "first")) 

exons = genes %>%
    separate_rows(exons, sep = ",") %>% 
    separate(exons, sep = "-", into = c("start", "end")) %>% 
    mutate(start = as.integer(start),
           end = as.integer(end),
           tpm=as.integer(tpm)) 

maxx = max(cov$end) + 1000
minn = 0


ref = ggplot(ref_genes, aes(start, ntx)) +
    ggplot2::geom_segment(data=ref_genes[ref_genes$direction=="last",],
                          aes(start, ntx, xend=end, yend=ntx),
                          arrow=arrow(ends="last",
                                      length = unit(0.05, "inches"))) +
    ggplot2::geom_segment(data=ref_genes[ref_genes$direction=="first",],
                          aes(start, ntx, xend=end, yend=ntx),
                          color = "black",
                          arrow=arrow(ends="first",
                                      length = unit(0.05, "inches"))) +
    geom_text(data=ref_genes, aes(end, ntx, label=name), size=3, hjust = 0) +
    theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.title.x = element_blank(),
          panel.background = element_rect(fill = "white"),
          plot.margin = margin( t = 0, b = 0,  unit = "pt")) +
    ylab("reference") +
    xlim(minn, maxx)

gcov  = ggplot(cov, aes(start, coverage)) +
    geom_line() +
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
          axis.title.x = element_blank(),
          panel.background = element_rect(fill = "white"),
          plot.margin = margin( b = 0,  unit = "pt")
    ) +
    xlim(minn, maxx)

stringtie = ggplot(exons, aes(start, ntx)) +
    ggplot2::geom_segment(aes(xend=end, yend=ntx,
                              color=log2(as.numeric(tpm))), size = 5) +
    ggplot2::geom_segment(data=genes[genes$direction=="last",],
                          aes(start, ntx, xend=end, yend=ntx),
                          arrow=arrow(ends="last",
                                      length = unit(0.05, "inches"))) +
    ggplot2::geom_segment(data=genes[genes$direction=="first",],
                          aes(start, ntx, xend=end, yend=ntx),
                          color = "black",
                          arrow=arrow(ends="first",
                                      length = unit(0.05, "inches"))) +
    geom_text(data=genes, aes(end, ntx, label=name), size=3, hjust = 0) +
    scale_color_continuous_tableau("log2(TPM)", palette = "Blue-Teal")  +
    theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(),
          panel.background = element_rect(fill = "white"),
          plot.margin = margin( t = 0, b = 0,  unit = "pt")) +
    ylab("predicted") +
    xlim(minn, maxx)

plot_grid(gcov, ref, stringtie, rel_heights = c(0.5,0.1,0.5), ncol=1, align = "v", axis = "lr" ) +
    ggsave2(output, height = 12, width = 9)

