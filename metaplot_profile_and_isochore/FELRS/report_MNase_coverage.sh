#!/bin/bash

## Generate the PPTX reports of the Felrs analysis
script_dir='/home/sebastien/analyses_SEBASTIEN/analyse_RNA-Seq_splicingFactors_2018-07/src/pdf'
wdir='/home/sebastien/analyses_SEBASTIEN/analyse_RNA-Seq_splicingFactors_2018-07/FELRS/GC-AT/res'
measure_window=50,500 #50,200
center_window=1000

auboeuf_vec='MNase_Seq,MNase_Seq_B1s,MNase_Seq_B2s,MNase_Seq_B3s'
good_public_vec='public_Raji_B3s_GSE36979,public_Raji_B4s_GSE36979,public_Raji_B5s_GSE36979,public_IMR90_B1s_GSE44985,public_WA09_hESC_B1s_GSE46461,public_WA09_INM_B1s_GSE46461,public_WA09_SMC_B1s_GSE46461,public_HUVEC_B1s_GSE53343,public_HUVEC_B2s_GSE53343,public_hESC_B1s_GSE59062,public_hESC_B2s_GSE59062,public_hESC_B3s_GSE59062,public_hFib_B1s_GSE59062,public_hFib_B2s_GSE59062,public_hFib_B3s_GSE59062,public_hFib_B4s_GSE59062,public_hiPSC_B1s_GSE59062,public_hiPSC_B2s_GSE59062,public_HCT116_B1s_GSE89871'
good_public_ChIP_H3='public_IMR90-prolif_H3_B1s_GSE36616,public_IMR90-senesc_H3_B1s_GSE36616,public_colon-carci_H3_B1s_GSE47190,public_LoVoDTT0_H3_B1s_GSE51290,public_MDA-MB-468_H3_B1s_GSE59176,public_HeLa_H3_B1s_GSE60526,public_A549_H3_B1s_GSE69885,public_AML-Gfi-36S_H3_B1s_GSE71254,public_AML-Gfi-36S_H3_B2s_GSE71254,public_Granulocyte_H3_B1s_GSE71809,public_K562_H3_B1s_GSE71809,public_IMR90_H3_B1s_GSE78141,public_LHCN-M2_H3_B1s_GSE78158,public_HL-60-S4_H3_B1s_GSE90992'

encode_chip_POLR2A='ENCSR000BGD_ENCFF593AJA_POLR2A_2,ENCSR000BGD_ENCFF644WID_POLR2A_1,ENCSR000BGO_ENCFF016XGN_POLR2A_1,ENCSR000BGO_ENCFF020STC_POLR2A_2,ENCSR000BIK_ENCFF348NIY_POLR2A_1,ENCSR000BIK_ENCFF734TAE_POLR2A_2,ENCSR000BKI_ENCFF095TCV_POLR2A_2,ENCSR000BKI_ENCFF759HDA_POLR2A_1,ENCSR000BMR_ENCFF166JUY_POLR2A_2,ENCSR000BMR_ENCFF749YKR_POLR2A_1,ENCSR000BQB_ENCFF423ZWY_POLR2A_2,ENCSR000BQB_ENCFF610EUL_POLR2A_1,ENCSR000DPA_ENCFF790YFZ_POLR2A_2,ENCSR000DPA_ENCFF971YJB_POLR2A_1,ENCSR025ZZA_ENCFF679ULL_POLR2A_1,ENCSR031KWR_ENCFF659PWI_POLR2A_1,ENCSR031TFS_ENCFF042PRZ_POLR2A_1,ENCSR031TFS_ENCFF345HKX_POLR2A_2,ENCSR033NHF_ENCFF845WMM_POLR2A_1,ENCSR038YKU_ENCFF436UEZ_POLR2A_1,ENCSR055WBT_ENCFF576AFZ_POLR2A_1,ENCSR068WNI_ENCFF418MEF_POLR2A_1,ENCSR091CSG_ENCFF551QFV_POLR2A_1,ENCSR103UPR_ENCFF172EHB_POLR2A_1,ENCSR107EUS_ENCFF202YYG_POLR2A_1,ENCSR122EYE_ENCFF758ZVL_POLR2A_1,ENCSR122LGV_ENCFF090VBQ_POLR2A_1,ENCSR132XRW_ENCFF391HIC_POLR2A_1,ENCSR142SQX_ENCFF277CSB_POLR2A_1,ENCSR147PYL_ENCFF606CVJ_POLR2A_1,ENCSR154GUK_ENCFF862FKS_POLR2A_2,ENCSR156UQO_ENCFF418EYD_POLR2A_1,ENCSR266XFE_ENCFF865RPA_POLR2A_1,ENCSR289VTP_ENCFF598SNZ_POLR2A_1,ENCSR299CAV_ENCFF728WWJ_POLR2A_1,ENCSR304IVU_ENCFF839CVW_POLR2A_1,ENCSR322JEO_ENCFF614PZE_POLR2A_1,ENCSR336YRS_ENCFF546CXV_POLR2A_1,ENCSR367UUC_ENCFF851WEB_POLR2A_1,ENCSR388QZF_ENCFF555CHY_POLR2A_1,ENCSR388QZF_ENCFF910GGT_POLR2A_2,ENCSR401ORC_ENCFF495LKH_POLR2A_1,ENCSR412QGD_ENCFF083XYK_POLR2A_1,ENCSR414SWG_ENCFF026EAA_POLR2A_1,ENCSR414TPL_ENCFF122VYV_POLR2A_1,ENCSR415MOW_ENCFF380ZRN_POLR2A_1,ENCSR431EHE_ENCFF539GCZ_POLR2A_1,ENCSR432GAO_ENCFF562HIW_POLR2A_1,ENCSR434MWV_ENCFF158TWJ_POLR2A_1,ENCSR472VBD_ENCFF953USH_POLR2A_1,ENCSR485BEB_ENCFF739LQP_POLR2A_1,ENCSR517FVL_ENCFF604WTQ_POLR2A_1,ENCSR543DUC_ENCFF864AQK_POLR2A_1,ENCSR576PII_ENCFF192SIL_POLR2A_1,ENCSR610EFT_ENCFF501FUQ_POLR2A_1,ENCSR667UWT_ENCFF488KWK_POLR2A_1,ENCSR674IEI_ENCFF791NBS_POLR2A_1,ENCSR691CPM_ENCFF366PVU_POLR2A_1,ENCSR696HTM_ENCFF757MYM_POLR2A_2,ENCSR699ZGH_ENCFF766ZCY_POLR2A_1,ENCSR724FCJ_ENCFF796DYJ_POLR2A_1,ENCSR729GXE_ENCFF034SBQ_POLR2A_1,ENCSR734ZOZ_ENCFF576LRL_POLR2A_1,ENCSR754DWU_ENCFF066UFH_POLR2A_1,ENCSR775VUA_ENCFF275QQZ_POLR2A_1,ENCSR803FAP_ENCFF666TUC_POLR2A_1,ENCSR825RBI_ENCFF284MDF_POLR2A_1,ENCSR877OYD_ENCFF158JDH_POLR2A_1,ENCSR888HXN_ENCFF574ZJP_POLR2A_1,ENCSR888OWU_ENCFF398DXB_POLR2A_1,ENCSR911JAX_ENCFF587CLH_POLR2A_1,ENCSR922GUA_ENCFF826GPP_POLR2A_1,ENCSR935XOT_ENCFF115NGL_POLR2A_1,ENCSR974HQI_ENCFF261ZGD_POLR2A_1,ENCSR979DMN_ENCFF755SZN_POLR2A_1,ENCSR995QNB_ENCFF874VTY_POLR2A_1,ENCSR999QVR_ENCFF913WSO_POLR2A_1'
encode_chip_POLR2B='ENCSR325RLL_ENCFF452WGO_POLR2B_1,ENCSR325RLL_ENCFF545XLF_POLR2B_2'
encode_chip_POLR2G='ENCSR283ZRI_ENCFF015NSS_POLR2G_1,ENCSR283ZRI_ENCFF665FQH_POLR2G_2,ENCSR854JES_ENCFF193HJC_POLR2G_1,ENCSR854JES_ENCFF738NEB_POLR2G_2'
encode_chip_eGFP_POLR2H='ENCSR400FSM_ENCFF003VRX_eGFP-POLR2H_1,ENCSR400FSM_ENCFF559WEI_eGFP-POLR2H_2'
encode_chip_POLR2AphosphoS2='ENCSR000DYF_ENCFF695AKF_POLR2AphosphoS2_2,ENCSR000DYF_ENCFF874ONS_POLR2AphosphoS2_1,ENCSR000DZK_ENCFF521QIX_POLR2AphosphoS2_1,ENCSR000DZK_ENCFF940CQB_POLR2AphosphoS2_2,ENCSR000EDX_ENCFF673FVQ_POLR2AphosphoS2_1,ENCSR000EDX_ENCFF863TEE_POLR2AphosphoS2_2,ENCSR000EGF_ENCFF528AUR_POLR2AphosphoS2_2,ENCSR000EGF_ENCFF987CNI_POLR2AphosphoS2_1'
encode_chip_POLR2AphosphoS5='ENCSR000BIA_ENCFF366RSZ_POLR2AphosphoS5_1,ENCSR000BIA_ENCFF647NER_POLR2AphosphoS5_2,ENCSR000BIC_ENCFF275IIJ_POLR2AphosphoS5_2,ENCSR000BIC_ENCFF786QHR_POLR2AphosphoS5_1,ENCSR000BIF_ENCFF327IZI_POLR2AphosphoS5_2,ENCSR000BIF_ENCFF769IAK_POLR2AphosphoS5_1,ENCSR000BIL_ENCFF105AYM_POLR2AphosphoS5_2,ENCSR000BIL_ENCFF914IHS_POLR2AphosphoS5_1,ENCSR000BML_ENCFF197KTW_POLR2AphosphoS5_1,ENCSR000BML_ENCFF718HZQ_POLR2AphosphoS5_2,ENCSR000BOV_ENCFF155ALO_POLR2AphosphoS5_2,ENCSR000BOV_ENCFF723FQX_POLR2AphosphoS5_1,ENCSR000BPA_ENCFF091UNG_POLR2AphosphoS5_2,ENCSR000BPA_ENCFF167CJH_POLR2AphosphoS5_1,ENCSR000BPC_ENCFF007MAV_POLR2AphosphoS5_1,ENCSR000BPC_ENCFF995HPP_POLR2AphosphoS5_2,ENCSR000BPI_ENCFF071KXB_POLR2AphosphoS5_2,ENCSR000BPI_ENCFF846ZDQ_POLR2AphosphoS5_1,ENCSR000BPL_ENCFF436YJR_POLR2AphosphoS5_2,ENCSR000BPL_ENCFF959LOI_POLR2AphosphoS5_1,ENCSR000BTL_ENCFF074QRH_POLR2AphosphoS5_1,ENCSR000BTL_ENCFF863IPV_POLR2AphosphoS5_2,ENCSR000BTW_ENCFF187WWV_POLR2AphosphoS5_1,ENCSR000BTW_ENCFF928JCP_POLR2AphosphoS5_2,ENCSR033FDW_ENCFF142CQJ_POLR2AphosphoS5_1,ENCSR170NMC_ENCFF760FNL_POLR2AphosphoS5_1,ENCSR350PUV_ENCFF946SLI_POLR2AphosphoS5_1,ENCSR369NGL_ENCFF433CWQ_POLR2AphosphoS5_2,ENCSR369NGL_ENCFF964CLX_POLR2AphosphoS5_1,ENCSR400WEK_ENCFF844DMB_POLR2AphosphoS5_1,ENCSR403USE_ENCFF502WCH_POLR2AphosphoS5_1,ENCSR404VWY_ENCFF753OGK_POLR2AphosphoS5_1,ENCSR442CIF_ENCFF204IHN_POLR2AphosphoS5_1,ENCSR442ZTI_ENCFF682YKY_POLR2AphosphoS5_1,ENCSR464CSO_ENCFF128YAW_POLR2AphosphoS5_1,ENCSR491PTJ_ENCFF079WAL_POLR2AphosphoS5_1,ENCSR633OEO_ENCFF691VMG_POLR2AphosphoS5_1,ENCSR645JVW_ENCFF645PDZ_POLR2AphosphoS5_1,ENCSR714EQS_ENCFF546RER_POLR2AphosphoS5_1,ENCSR714JRB_ENCFF460SOG_POLR2AphosphoS5_1,ENCSR844PVS_ENCFF986GSI_POLR2AphosphoS5_1,ENCSR861XGM_ENCFF215PWS_POLR2AphosphoS5_1,ENCSR889GGV_ENCFF285UTG_POLR2AphosphoS5_2,ENCSR889GGV_ENCFF701CDH_POLR2AphosphoS5_1,ENCSR893MYW_ENCFF867HCH_POLR2AphosphoS5_1,ENCSR916JAC_ENCFF642XTK_POLR2AphosphoS5_1,ENCSR919UCY_ENCFF216SHO_POLR2AphosphoS5_1,ENCSR920FID_ENCFF194WKG_POLR2AphosphoS5_1,ENCSR930CPA_ENCFF547WRK_POLR2AphosphoS5_1,ENCSR936JHB_ENCFF619RKX_POLR2AphosphoS5_1,ENCSR939FGB_ENCFF506QKU_POLR2AphosphoS5_1,ENCSR950CUQ_ENCFF245MLE_POLR2AphosphoS5_1,ENCSR959UQA_ENCFF563XUY_POLR2AphosphoS5_1,ENCSR978LQC_ENCFF196PEM_POLR2AphosphoS5_1,ENCSR986QAR_ENCFF848JVN_POLR2AphosphoS5_1'


## specific variable setters
Hela_fc () {
  ref_list='Hela'
  pptx_suffix=GCbias_Nbs_${ref_list}
  measure_list='public_HeLa_B1s_GSE65644,public_HeLaS3_B1s_GSE14022,public_HeLaS3_B1s_GSE42951,public_HeLaS3_B2s_GSE42951',${auboeuf_vec},${good_public_vec}
  #
}

K562_fc () {
  ref_list='K562'
  pptx_suffix=GCbias_Nbs_${ref_list}
  measure_list='public_K562_5U_B1s_GSE78984,public_K562_20U_B1s_GSE78984,public_K562_80U_B1s_GSE78984,public_K562_300U_B1s_GSE78984,public_K562_B1s_GSE70920,public_K562_low_B1s_GSE84474,public_K562_medium_B1s_GSE84474,public_K562_high_B1s_GSE84474',${auboeuf_vec},${good_public_vec}
}

all_fc () {
  ref_list='GC-AT'
  pptx_suffix=GCbias_Nbs_${ref_list}
  measure_list='public_HeLa_B1s_GSE65644,public_HeLaS3_B1s_GSE14022,public_HeLaS3_B1s_GSE42951,public_HeLaS3_B2s_GSE42951,public_K562_5U_B1s_GSE78984,public_K562_20U_B1s_GSE78984,public_K562_80U_B1s_GSE78984,public_K562_300U_B1s_GSE78984,public_K562_B1s_GSE70920,public_K562_low_B1s_GSE84474,public_K562_medium_B1s_GSE84474,public_K562_high_B1s_GSE84474',${auboeuf_vec},${good_public_vec}
}

good_MNase-Seq_fc () {
  ref_list='GC-AT'
  pptx_suffix=GCbias_Nbs_${ref_list}
  measure_list=${auboeuf_vec},${good_public_vec}
}

CCE_good_MNase-Seq_fc () {
  ref_list='CCE'
  pptx_suffix=CCE_Nbs_${ref_list}
  measure_list=${auboeuf_vec},${good_public_vec}
}

good_ChIP_H3_fc () {
  ref_list='GC-AT'
  pptx_suffix=GCbias_Nbs_${ref_list}_H3
  measure_list=${good_public_ChIP_H3}
}

good_ChIP_POLR2A_fc () {
  ref_list='GC-AT'
  pptx_suffix=GCbias_Nbs_${ref_list}_POLR2A
  measure_list=${encode_chip_POLR2A}
}

good_ChIP_POLR2B_fc () {
  ref_list='GC-AT'
  pptx_suffix=GCbias_Nbs_${ref_list}_POLR2B
  measure_list=${encode_chip_POLR2B}
}

good_ChIP_POLR2G_fc () {
  ref_list='GC-AT'
  pptx_suffix=GCbias_Nbs_${ref_list}_POLR2G
  measure_list=${encode_chip_POLR2A}
}

good_ChIP_eGFP-POLR2H_fc () {
  ref_list='GC-AT'
  pptx_suffix=GCbias_Nbs_${ref_list}_eGFP-POLR2H
  measure_list=${encode_chip_POLR2A}
}

good_ChIP_POLR2AphosphoS2_fc () {
  ref_list='GC-AT'
  pptx_suffix=GCbias_Nbs_${ref_list}_POLR2AphosphoS2
  measure_list=${encode_chip_POLR2A}
}

good_ChIP_POLR2AphosphoS5_fc () {
  ref_list='GC-AT'
  pptx_suffix=GCbias_Nbs_${ref_list}_POLR2AphosphoS5
  measure_list=${encode_chip_POLR2A}
}

CCE_good_ChIP_H3_fc () {
  ref_list='CCE'
  pptx_suffix=CCE_Nbs_${ref_list}_H3
  measure_list=${good_public_ChIP_H3}
}

GC-AT_BF_fc () {
  ref_list='GC-AT'
  pptx_suffix=GCbias_Nbs_${ref_list}
  measure_list='A,T,G,C,GC,AT,AG,TC,AC,TG'
}

CCE_BF_fc () {
  ref_list='CCE'
  pptx_suffix=CCE_Nbs_${ref_list}
  measure_list='A,T,G,C,GC,AT,AG,TC,AC,TG'
}


#### MAIN

## on Exons
# spec_var_set_list='Hela_fc K562_fc'
#
# for extra_suf in '' '_byExp' '_PPT_bias_byExp'; do
#   for spec_var_set in ${spec_var_set_list}; do
#     ${spec_var_set}
#
#     python3 \
#     ${script_dir}/pptx_measure_multiAnal_builder.py \
#     ${pptx_suffix}${extra_suf} \
#     ${wdir}/${ref_list}_nbup${extra_suf},${wdir}/${ref_list}${extra_suf},${wdir}/${ref_list}_nbdown${extra_suf} \
#     ${measure_window} \
#     ${measure_list} \
#     coverage \
#     --out ${wdir}/../pptx_reports \
#     --res-name n-1,n,n+1 \
#     ;
#   done
# done
#
#
##
all_fc
python3 \
${script_dir}/pptx_measure_multiAnal_builder.py \
${pptx_suffix} \
${wdir}/${ref_list}_nbup,${wdir}/${ref_list},${wdir}/${ref_list}_nbdown \
${measure_window} \
${measure_list} \
coverage \
--out ${wdir}/../pptx_reports \
--res-name n-1,n,n+1 \
;


good_MNase-Seq_fc
python3 \
${script_dir}/pptx_measure_multiAnal_builder.py \
${pptx_suffix} \
${wdir}/${ref_list}_nbup,${wdir}/${ref_list},${wdir}/${ref_list}_nbdown \
${measure_window} \
${measure_list} \
coverage \
--out ${wdir}/../pptx_reports \
--res-name n-1,n,n+1 \
;

good_MNase-Seq_fc
python3 \
${script_dir}/pptx_measure_multiAnal_builder.py \
${pptx_suffix}_center \
${wdir}/${ref_list}_nbup,${wdir}/${ref_list},${wdir}/${ref_list}_nbdown \
${center_window} \
${measure_list} \
coverage \
--out ${wdir}/../pptx_reports \
--res-name n-1,n,n+1 \
--center \
;

## Pol2
for target in "eGFP_POLR2H" "POLR2A" "POLR2AphosphoS2" "POLR2AphosphoS5" "POLR2B" "POLR2G"; do
  good_ChIP_${target}_fc
  python3 \
  ${script_dir}/pptx_measure_multiAnal_builder.py \
  ${pptx_suffix} \
  ${wdir}/${ref_list}_nbup,${wdir}/${ref_list},${wdir}/${ref_list}_nbdown \
  ${measure_window} \
  ${measure_list} \
  coverage \
  --out ${wdir}/../pptx_reports \
  --res-name n-1,n,n+1 \
  ;

  python3 \
  ${script_dir}/pptx_measure_multiAnal_builder.py \
  ${pptx_suffix}_center \
  ${wdir}/${ref_list}_nbup,${wdir}/${ref_list},${wdir}/${ref_list}_nbdown \
  ${center_window} \
  ${measure_list} \
  coverage \
  --out ${wdir}/../pptx_reports \
  --res-name n-1,n,n+1 \
  --center \
  ;
done

CCE_good_MNase-Seq_fc
python3 \
${script_dir}/pptx_measure_multiAnal_builder.py \
${pptx_suffix}_center \
${wdir}/${ref_list} \
${center_window} \
${measure_list} \
coverage \
--out ${wdir}/../pptx_reports \
--res-name n \
--center \
;

CCE_good_MNase-Seq_fc
python3 \
${script_dir}/pptx_measure_multiAnal_builder.py \
${pptx_suffix} \
${wdir}/${ref_list} \
${measure_window} \
${measure_list} \
coverage \
--out ${wdir}/../pptx_reports \
--res-name n \
;


## ChIP H3
good_ChIP_H3_fc
python3 \
${script_dir}/pptx_measure_multiAnal_builder.py \
${pptx_suffix} \
${wdir}/${ref_list} \
${measure_window} \
${measure_list} \
coverage \
--out ${wdir}/../pptx_reports \
--res-name n \
;

CCE_good_ChIP_H3_fc
python3 \
${script_dir}/pptx_measure_multiAnal_builder.py \
${pptx_suffix}_center \
${wdir}/${ref_list} \
${center_window} \
${measure_list} \
coverage \
--out ${wdir}/../pptx_reports \
--res-name n \
--center \
;

CCE_good_ChIP_H3_fc
python3 \
${script_dir}/pptx_measure_multiAnal_builder.py \
${pptx_suffix} \
${wdir}/${ref_list} \
${measure_window} \
${measure_list} \
coverage \
--out ${wdir}/../pptx_reports \
--res-name n \
;


## Base Fraction
for ext_win in 1 20;do
  GC-AT_BF_fc
  python3 \
  ${script_dir}/pptx_measure_multiAnal_builder.py \
  ${pptx_suffix}_w${ext_win} \
  ${wdir}/${ref_list}_nbup,${wdir}/${ref_list},${wdir}/${ref_list}_nbdown \
  ${measure_window},${ext_win} \
  ${measure_list} \
  base_fraction \
  --out ${wdir}/../pptx_reports \
  --res-name n-1,n,n+1 \
  ;
done

for ext_win in 1 20;do
  GC-AT_BF_fc
  python3 \
  ${script_dir}/pptx_measure_multiAnal_builder.py \
  ${pptx_suffix}_w${ext_win} \
  ${wdir}/${ref_list} \
  ${measure_window},${ext_win} \
  ${measure_list} \
  base_fraction \
  --out ${wdir}/../pptx_reports \
  --res-name n \
  ;
done


###############
## on Nucleosomes
wdir='/home/sebastien/analyses_SEBASTIEN/analyse_RNA-Seq_splicingFactors_2018-07/FELRS/nucleosome/res'

nuc_fc () {
  ref_list='Nuc'
  pptx_suffix=Nuc_${ref_list}
  measure_list=${auboeuf_vec}
}

nuc_h3_fc () {
  ref_list='Nuc'
  pptx_suffix=Nuc_${ref_list}_H3
  measure_list=${good_public_ChIP_H3}
}

#### MAIN

nuc_fc
python3 \
${script_dir}/pptx_measure_multiAnal_builder.py \
${pptx_suffix}_center \
${wdir}/${ref_list} \
${center_window} \
${measure_list} \
coverage \
--out ${wdir}/../pptx_reports \
--res-name n \
--center \
;

nuc_h3_fc
python3 \
${script_dir}/pptx_measure_multiAnal_builder.py \
${pptx_suffix}_center \
${wdir}/${ref_list} \
${center_window} \
${measure_list} \
coverage \
--out ${wdir}/../pptx_reports \
--res-name n \
--center \
;


####
