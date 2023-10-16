# synteny_plotter

A tool to visualise synteny relationships between genomes.

Testing run

```
Rscript scripts/generate_synteny_plot.R -busco1 test_data/Melitaea_cinxia.tsv -busco2 test_data/Vanessa_cardui.tsv -chrom1 test_data/Melitaea_cinxia_info.tsv -chrom2 test_data/Vanessa_cardui_info.tsv -o test_data/test
```

## automated generator of chromoamal order

The following reads testing files, produces order of chromosomes for Vanessa (`test_data/Vanessa_cardui_info_generated.tsv`)

```
Rscript scripts/dev_generate_chromosome_file.R
```

Now I would like to compare the manual and generated versions

```
Rscript scripts/generate_synteny_plot.R  -busco1 test_data/Melitaea_cinxia.tsv -busco2 test_data/Vanessa_cardui.tsv -chrom1 test_data/Melitaea_cinxia_info.tsv -chrom2 test_data/Vanessa_cardui_curated_info.tsv -o test_data/test_curated

Rscript scripts/generate_synteny_plot.R -busco1 test_data/Melitaea_cinxia.tsv -busco2 test_data/Vanessa_cardui.tsv -chrom1 test_data/Melitaea_cinxia_info.tsv -chrom2 test_data/Vanessa_cardui_info_generated.tsv -o test_data/test_generated
```

The order seems to work well.

