import allel

vcf_file = "/mnt/storage13/ahri/training_FA5/admixture/012026_pf_mgen_ug_india_eritrea_ethiopia_popgen.filt.csq.bi.GT.miss0.2.vqslod.filt.snps.erit.vcf.gz"
zarr_output = "/mnt/storage13/ahri/training_FA5/admixture/012026_pf_mgen_ug_india_eritrea_ethiopia_popgen.filt.csq.bi.GT.miss0.2.vqslod.filt.snps.erit.zarr"

allel.vcf_to_zarr(
    vcf_file,
    zarr_output,
    fields='*',
    overwrite=True
)