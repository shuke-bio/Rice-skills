#!/usr/bin/env python3
###############################################################################
# 04_scvi_cellchat.py -- scVI integration + CellChat plan/runner
#
# 说明：CellChat 为 R 包，脚本生成 R 计划；scVI 若已安装且 --run 则执行。
#
# 用法：
#   python scripts/single_cell/04_scvi_cellchat.py --h5ad data.h5ad --batch-key batch \
#       --celltype-key cell_type --outdir results/single_cell/04_scvi --run
#   python scripts/single_cell/04_scvi_cellchat.py --demo --outdir results/single_cell/04_scvi_demo
###############################################################################
import argparse, os, sys, datetime

def write_plans(outdir, h5ad, batch, celltype):
    py = f'''# scVI integration plan
import scanpy as sc, scvi, anndata as ad
adata = ad.read_h5ad("{h5ad}")
scvi.model.SCVI.setup_anndata(adata, batch_key="{batch or 'batch'}")
model = scvi.model.SCVI(adata, n_latent=30)
model.train(max_epochs=200)
adata.obsm["X_scVI"] = model.get_latent_representation()
sc.pp.neighbors(adata, use_rep="X_scVI")
sc.tl.leiden(adata, resolution=0.8, key_added="scvi_cluster")
adata.obs[["scvi_cluster"]].to_csv("{outdir}/scvi_clusters.csv")
'''
    open(os.path.join(outdir, "scvi_plan.py"), "w", encoding="utf-8").write(py)
    r = f'''# CellChat plan (R)
library(CellChat)
data.input <- read.csv("{outdir}/cellchat_expression.csv", row.names=1)
meta <- read.csv("{outdir}/cellchat_metadata.csv", row.names=1)
cellchat <- createCellChat(object=data.input, meta=meta, group.by="{celltype or 'cell_type'}")
# cellchat@DB <- CellChatDB.human  # 请替换为植物/自定义 LR 库
cellchat <- subsetData(cellchat)
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- computeCommunProb(cellchat)
cellchat <- filterCommunication(cellchat, min.cells=10)
df.net <- subsetCommunication(cellchat)
write.csv(df.net, "{outdir}/cellchat_communications.csv", row.names=FALSE)
'''
    open(os.path.join(outdir, "cellchat_plan.R"), "w", encoding="utf-8").write(r)

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--h5ad"); ap.add_argument("--batch-key", default=None)
    ap.add_argument("--celltype-key", default=None)
    ap.add_argument("--outdir", default="results/single_cell/04_scvi")
    ap.add_argument("--run", action="store_true"); ap.add_argument("--demo", action="store_true")
    a = ap.parse_args()
    os.makedirs(a.outdir, exist_ok=True)
    if a.demo:
        a.h5ad = "demo.h5ad"
    if not a.h5ad: sys.exit("Need --h5ad or --demo")
    write_plans(a.outdir, a.h5ad, a.batch_key, a.celltype_key)
    executed = False
    if a.run and not a.demo:
        try:
            import anndata as ad, scanpy as sc, scvi
            adata = ad.read_h5ad(a.h5ad)
            scvi.model.SCVI.setup_anndata(adata, batch_key=a.batch_key or None)
            model = scvi.model.SCVI(adata, n_latent=30)
            model.train(max_epochs=200)
            adata.obsm["X_scVI"] = model.get_latent_representation()
            sc.pp.neighbors(adata, use_rep="X_scVI")
            sc.tl.leiden(adata, resolution=0.8, key_added="scvi_cluster")
            adata.obs[["scvi_cluster"]].to_csv(os.path.join(a.outdir, "scvi_clusters.csv"))
            executed = True
        except Exception as e:
            print(f"[scvi] not executed: {e}")
    open(os.path.join(a.outdir, "run_info.txt"), "w", encoding="utf-8").write(
        f"date={datetime.datetime.now().isoformat(timespec='seconds')}\n"
        f"h5ad={a.h5ad}\nexecuted={executed}\ncellchat=plan only (R)\n")
    print(f"[scvi] plans -> {a.outdir}/scvi_plan.py, cellchat_plan.R (executed={executed})")

if __name__ == "__main__":
    main()
