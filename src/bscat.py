"""
bscat.py
Scans in Boundary Condition to generate wavefunctions with various energies
"""
#-----------------------------------------------------------------------
import numpy as np
import json
#-----------------------------------------------------------------------
from utility import utility_t
from deck import deck_t,GenerateOptFile
from wavefunction import wavefunction_t
#-----------------------------------------------------------------------
"""
---
Single Channel Logic
---
1. Assume that core/target is optimized
2. Move A-1 parameters into A body input deck
3. Optimize with inital bscat
4. Initialize scan log
5. Optimize with +- initial_delta_bscat
6. Update log with +-
7. Compute db/de (central finite difference)
8. Scan Increasing bscat from initial
9. Scan Decreasing bscat from initial
---
Scan Logic
---
1. direction = +-
2. bscat_next = bscat_inital + direction * |db_de| * de
3. optimize with bscat_next
4. update log
5. update db_de (foward/backward finite difference)
6. check bounds, emin,emax,db_de_max
7. todo: if db_de > db_de_max and e_current < emax => node_count+=1 bscat = 2.5
   ***care needed in wse search
---
Optimize Logic
---
1. INPUT: bscat + deck + spatial symmetry index, logname
2. Scan WSE for the Spatial Symmetry (defualt: delta=0.5,num_walks=5,samples_per_walk=5)
3. Optimize correlations that only affect pairs/triples that include scattering nucleon (without wse)
4. Optimize correlations that only affect pairs/triples that include scattering nucleon (with wse)
"""
def SingleChannelOptimize(bscat:float, \
                          work_dir:str, \
                          label:str, \
                          scatter:wavefunction_t,\
                          ssi:int,\
                          opt_scale:float,\
                          ecore:float,\
                          vcore:float):
#-----------------------------------------------------------------------
    BREAK="="*72
#-----------------------------------------------------------------------
    opt_wse_file_name=f"\'{work_dir}opt/wse.opt\'"
    wse_i=[{"ss":True,"ss_idx":ssi,"key":"WSE","scale":0, "flat":0.5}]
#-----------------------------------------------------------------------
    instructions=[{"ss":False,"key":"QSSP1","scale":opt_scale, "flat":0.0}\
                 ,{"ss":False,"key":"QSSP2","scale":opt_scale, "flat":0.0}\
                 ,{"ss":True,"ss_idx":ssi,"key":"SPU"  ,"scale":opt_scale, "flat":0.0}\
                 ,{"ss":True,"ss_idx":ssi,"key":"SPV"  ,"scale":opt_scale, "flat":0.0}\
                 ,{"ss":True,"ss_idx":ssi,"key":"SPA"  ,"scale":opt_scale, "flat":0.0}\
                 ,{"ss":True,"ss_idx":ssi,"key":"SPB"  ,"scale":opt_scale, "flat":0.0}\
                 ,{"ss":True,"ss_idx":ssi,"key":"SPC"  ,"scale":opt_scale, "flat":0.0}\
                 ,{"ss":True,"ss_idx":ssi,"key":"SPK"  ,"scale":opt_scale, "flat":0.0}\
                 ,{"ss":True,"ss_idx":ssi,"key":"SPL"  ,"scale":opt_scale, "flat":0.0}\
                 ,{"ss":True,"ss_idx":ssi,"key":"WSR"  ,"scale":opt_scale, "flat":0.0}\
                 ,{"ss":True,"ss_idx":ssi,"key":"WSA"  ,"scale":opt_scale, "flat":0.0}]
#-----------------------------------------------------------------------
    instructions_all=instructions+wse_i
#-----------------------------------------------------------------------
    opt_wse=GenerateOptFile(scatter.PARAMS, scatter.DK, opt_wse_file_name,wse_i)
#-----------------------------------------------------------------------
    num_samples=scatter.CTRL.NUM_OPT_EVALUATIONS
#-----------------------------------------------------------------------
    scatter.DK=deck_t(scatter.PARAMS,scatter.DK.FILE_NAME)
    scatter.DK.FILE_NAME=f"\'{work_dir}temp.dk\'"
    scatter.DK.SS[ssi].BSCAT=str(bscat)
    scatter.DK.Write(scatter.PARAMS,scatter.DK.FILE_NAME)
    scatter.CTRL.INPUT_BRA.DECK_FILE=scatter.DK.FILE_NAME
#-----------------------------------------------------------------------
    dk_name=f"\'{work_dir}dk/{label}.{bscat:.4f}.dk\'"
#-----------------------------------------------------------------------
    print(f" 🔥 BSCAT = {bscat:.4f} 🔥")
    print(BREAK)
    log=f"{work_dir}logs/{label}.start.{bscat:.4f}"
    print(f"BEGIN INITIAL EVALUATION: {log}.energy")
    e,v=scatter.Evaluate(True,log)
    e_start=e-ecore
    v_start=np.sqrt(v**2+vcore**2)
    print(f" ⚛ EREL = {e_start:.4f} +- {v_start:.4f}")
#-----------------------------------------------------------------------
    scatter.CTRL.NUM_OPT_EVALUATIONS="5"
    scatter.CTRL.NUM_OPT_WALKS="5"
#-----------------------------------------------------------------------
    log=f"{work_dir}logs/{label}.wse.{bscat:.4f}"
    print(f"BEGIN WSE SCAN: {log}.optimize")
    e,v=scatter.Optimize(opt_wse,dk_name,True,log)
    wse_erel=e-ecore
    wse_vrel=np.sqrt(v**2+vcore**2)
    print(f"OPTIMAL WSE = {float(scatter.DK.SS[ssi].WSE):.4f}")
    print(f" ⚛ EREL = {wse_erel:.4f} +- {wse_vrel:.4f}")
    print(f"WSE OPTIMIZATION LOWERED ENERGY BY: {wse_erel-e_start:.4f} MeV")
#-----------------------------------------------------------------------
    scatter.CTRL.NUM_OPT_EVALUATIONS=num_samples
    scatter.CTRL.NUM_OPT_WALKS="1"
#-----------------------------------------------------------------------
    log=f"{work_dir}logs/{label}.corr.{bscat:.4f}"
    print(f"BEGIN CORRELATION OPTIMIZATION: {log}.optimize")
    opt_corr_file_name=f"\'{work_dir}opt/{label}.corr.{bscat:.4f}.opt\'"
    opt_corr=GenerateOptFile(scatter.PARAMS, scatter.DK, opt_corr_file_name,instructions)
    opt_corr.UpdateFloats(scatter.PARAMS,6)
    e,v=scatter.Optimize(opt_corr,dk_name,True,log)
    corr_erel=e-ecore
    corr_vrel=np.sqrt(v**2+vcore**2)
    print(f" ⚛ EREL = {corr_erel:.4f} +- {corr_vrel:.4f}")
    print(f"CORRELATION OPTIMIZATION LOWERED ENERGY BY: {corr_erel-wse_erel:.4f} MeV")
#-----------------------------------------------------------------------
    log=f"{work_dir}logs/{label}.all.{bscat:.4f}"
    print(f"BEGIN CORRELATION OPTIMIZATION: {log}.optimize")
    opt_all_file_name=f"\'{work_dir}opt/{label}.all.{bscat:.4f}.opt\'"
    opt_all=GenerateOptFile(scatter.PARAMS, scatter.DK, opt_all_file_name,instructions_all)
    opt_all.UpdateFloats(scatter.PARAMS,6)
    e,v=scatter.Optimize(opt_all,dk_name,True,log)
    all_erel=e-ecore
    all_vrel=np.sqrt(v**2+vcore**2)
    print(f" ⚛ EREL = {all_erel:.4f} +- {all_vrel:.4f}")
    print(f"ALL OPTIMIZATION LOWERED ENERGY BY: {all_erel-corr_erel:.4f} MeV")
    print(BREAK)
#-----------------------------------------------------------------------
    dat = {
            "TYPE": "VMC",
            "BSCAT": bscat,
            "MU": float(scatter.PARAMS.NPART-1)/float(scatter.PARAMS.NPART),
            "L": int(scatter.DK.SS[ssi].LQNUM),
            "R": float(scatter.CTRL.BOX_SIZE),
            "NODES": int(scatter.DK.SS[ssi].LNODES),
            "EREL": all_erel,
            "VREL": float(all_vrel),
            "WSE": float(scatter.DK.SS[ssi].WSE),
            "DECK_PATH": dk_name
          }
    print(dat)
    print(BREAK)
    return dat

def DirectionalScan(direction,\
                    dk_file_i,\
                    bscat_i,\
                    erel_i,\
                    db_de_i,\
                    util:utility_t,\
                    label:str,\
                    scatter:wavefunction_t,\
                    ssi:int,\
                    ecore:float,vcore:float):
#-----------------------------------------------------------------------
    BREAK="="*72
#-----------------------------------------------------------------------
    scatter.DK.FILE_NAME=dk_file_i
    bscat_prev=bscat_i
    erel_prev=erel_i
    db_de=db_de_i
    do_scan=True
    add_node=False
    de=0
    idx=0
    dscan=[]
#-----------------------------------------------------------------------
    while (do_scan): #direction=-1
#-----------------------------------------------------------------------
        db=db_de*util.DELTA_ENERGY
        print(f"db: {db:.4f} :: de: {de:.4f} :: db/de :: {direction*db_de:.4f}")
        print(BREAK)
        bscat=bscat_prev+db*direction
        dscan.append(SingleChannelOptimize(bscat,util.WORKING_DIR,label,scatter,ssi,util.OPT_SCALE,ecore,vcore))
#-----------------------------------------------------------------------
        de=dscan[idx]["EREL"]-erel_prev
#-----------------------------------------------------------------------
        db_de=db/abs(de)
#-----------------------------------------------------------------------
        bscat_prev=bscat
        erel_prev=dscan[idx]["EREL"]
#-----------------------------------------------------------------------
        idx+=1
#-----------------------------------------------------------------------
        out_of_bounds = (erel_prev < util.ENERGY_LOWER_BOUND) or (erel_prev > util.ENERGY_UPPER_BOUND)
        add_node = (de > 0 ) and (db_de > util.MAX_BSCAT_SLOPE) and (bscat*direction > 0.) #increasing energy and db_de above threshold
#-----------------------------------------------------------------------
        if out_of_bounds:
            do_scan = False
            print(" 🏁 ENERGY OUT OF BOUNDS")
            print(BREAK)
        if add_node: do_scan = False
#-----------------------------------------------------------------------
    if add_node and not out_of_bounds: 
        print(" 🫚 ADDING NODE TO PHI")
        print(BREAK)
        scatter.DK.SS[ssi].LNODES=str(int(scatter.DK.SS[ssi].LNODES)+1)
        scatter.DK.FILE_NAME=f"\'{util.WORKING_DIR}temp.dk\'"
        scatter.DK.Write(scatter.PARAMS,scatter.DK.FILE_NAME)
        for bscat in [10.,5.,2.5,1.25,0.625]:
            dscan.append(SingleChannelOptimize(bscat*(-1)*direction\
                ,util.WORKING_DIR,label,scatter,ssi,util.OPT_SCALE,ecore,vcore))
#-----------------------------------------------------------------------
    return dscan
#-----------------------------------------------------------------------
def SingleChannelScan(util:utility_t,\
                      label:str,\
                      scatter:wavefunction_t,\
                      ssi:int,\
                      ecore:float,vcore:float):
#-----------------------------------------------------------------------
    scan=[]
    BREAK="="*72
#-----------------------------------------------------------------------
# INITIAL BSCATS
#-----------------------------------------------------------------------
    bscat=util.INITIAL_BSCAT
    scan.append(SingleChannelOptimize(bscat,util.WORKING_DIR,label,scatter,ssi,util.OPT_SCALE,ecore,vcore))
    bscat=util.INITIAL_BSCAT+util.INITIAL_DELTA_BSCAT
    scan.append(SingleChannelOptimize(bscat,util.WORKING_DIR,label,scatter,ssi,util.OPT_SCALE,ecore,vcore))
    bscat=util.INITIAL_BSCAT-util.INITIAL_DELTA_BSCAT
    scatter.DK.FILE_NAME=scan[0]["DECK_PATH"]
    scan.append(SingleChannelOptimize(bscat,util.WORKING_DIR,label,scatter,ssi,util.OPT_SCALE,ecore,vcore))
#-----------------------------------------------------------------------
    db_de_i=(2.*util.INITIAL_DELTA_BSCAT)/abs(scan[2]["EREL"]-scan[1]["EREL"])
#-----------------------------------------------------------------------
    print("📈 INCREASING BSCAT SCAN")
    print(BREAK)
#-----------------------------------------------------------------------
    scan+=DirectionalScan(1,\
                          scan[1]["DECK_PATH"],\
                          util.INITIAL_BSCAT,\
                          scan[0]["EREL"],\
                          db_de_i,\
                          util,\
                          label,\
                          scatter,\
                          ssi,\
                          ecore,vcore)
#-----------------------------------------------------------------------
    print("📉 DECREASING BSCAT SCAT")
    print(BREAK)
#-----------------------------------------------------------------------
    scan+=DirectionalScan(-1,\
                          scan[2]["DECK_PATH"],\
                          util.INITIAL_BSCAT,\
                          scan[0]["EREL"],\
                          db_de_i,\
                          util,\
                          label,\
                          scatter,\
                          ssi,\
                          ecore,vcore)
#-----------------------------------------------------------------------
    scan_sorted = sorted(scan, key=lambda r: r["EREL"])
    with open(f"{util.WORKING_DIR}{label}.vmc.json", "w") as f:
        json.dump(scan_sorted, f, indent=2)
#-----------------------------------------------------------------------
    return 0