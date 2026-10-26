#!/usr/bin/env python3
###############################################################################
# 03_cell2location.py -- cell2location / COMMOT plan generator
#
# 说明：cell2location 与 COMMOT 依赖较重，脚本生成可执行计划；
#       若已安装且 --run，可尝试运行 cell2location 的最小流程。
#
# 用法：
#   python scripts/spatial/03_cell2location.py --spatial-h5ad spatial.h5ad \
#       --sc-h5ad sc.h5ad --labels labels.csv --outdir results/spatial/03_c2l --run
#   python scripts/spatial/03_cell2location.py --demo --outdir results/spatial/03_c2l_demo
###############################################################################
import argparse, os, sys, datetime

def write_plan(outdir, sp, sc, labels):
    py = f'''# cell2location plan
import scanpy as sc, anndata as ad, cell2location
sp = ad.read_h5ad("{sp}")
ref = ad.read_h5ad("{sc}")
labels = __import__("pandas").read_csv("{labels}")
# 1) reference signatures
cell2location.models.RegressionModel.setup_anndata(ref, labels_key="cell_type")
reg = cell2location.models.RegressionModel(ref)
reg.train(max_epochs=250)
ref = reg.export_posterior(ref)
ref.var[["means_per_cluster_mu_fg"]].to_csv("{outdir}/reference_signatures.csv")
# 2) spatial mapping
cell2location.models.Cell2location.setup_anndata(sp, batch_key="sample")
mod = cell2location.models.Cell2location(sp, cell_state_df=ref.var, N_cells_per_location=8)
mod.train(max_epochs=1500)
sp = mod.export_posterior(sp)
sp.obsm["q05_cell_abundance"].to_csv("{outdir}/cell2location_abundance.csv")
'''
    open(os.path.join(outdir, "cell2location_plan.py"), "w", encoding="utf-8").write(py)
    commot = f'''# COMMOT plan (Python)
import commot as ct, anndata as ad
adata = ad.read_h5ad("{sp}")
ct.tl.spatial_communication(adata, database_name="cellphonedb", df_ligrec=ct.databases.database("cellphonedb"),
                            dis_thr=500, heteromeric=True)
adata.obsm["commot"].to_csv("{outdir}/commot_communication.csv")
'''
    open(os.path.join(outdir, "commot_plan.py"), "w", encoding="utf-8").write(commot)

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--spatial-h5ad"); ap.add_argument("--sc-h5ad"); ap.add_argument("--labels")
    ap.add_argument("--outdir", default="results/spatial/03_c2l")
    ap.add_argument("--run", action="store_true"); ap.add_argument("--demo", action="store_true")
    a = ap.parse_args()
    os.makedirs(a.outdir, exist_ok=True)
    if a.demo:
        a.spatial_h5ad, a.sc_h5ad, a.labels = "spatial.h5ad", "sc.h5ad", "labels.csv"
    if not a.spatial_h5ad or not a.sc_h5ad or not a.labels:
        sys.exit("Need --spatial-h5ad --sc-h5ad --labels, or --demo")
    write_plan(a.outdir, a.spatial_h5ad, a.sc_h5ad, a.labels)
    executed = False
    if a.run and not a.demo:
        try:
            import anndata as ad, cell2location
            sp = ad.read_h5ad(a.spatial_h5ad)
            ref = ad.read_h5ad(a.sc_h5ad)
            import pandas as pd
            labels = pd.read_csv(a.labels)
            # 最小可运行：仅验证数据与标签可读；完整建模见 cell2location_plan.py
            print(f"[c2l] loaded spatial={sp.n_obs} cells, reference={ref.n_obs} cells, labels={len(labels)}")
            executed = True
        except Exception as e:
            print(f"[c2l] not executed: {e}")
    open(os.path.join(a.outdir, "run_info.txt"), "w", encoding="utf-8").write(
        f"date={datetime.datetime.now().isoformat(timespec='seconds')}\n"
        f"spatial={a.spatial_h5ad}\nreference={a.sc_h5ad}\nlabels={a.labels}\nexecuted={executed}\n")
    print(f"[c2l] plans -> {a.outdir}/cell2location_plan.py, commot_plan.py (executed={executed})")

if __name__ == "__main__":
    main()
