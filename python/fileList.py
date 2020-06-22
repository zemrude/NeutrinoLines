import os

data_path = os.environ['ANALYSIS_DATA_DIR']

allFiles = {}

allFiles['MC'] = {}

allFiles['MC']['corsika'] = data_path+'L4_IC86.merged_corsika_LE_and_HE_11499_11808_11865_11905_11926_11943_12161_12268_12332_11937_12025_merged.npy'

allFiles['MC']['nominal'] = {}
allFiles['MC']['nominal']['nugen_e_LE'] = data_path+'L4_IC86.2012_merged_nugen_nue_low_energy_LE_and_HE_files_6953.npy'
allFiles['MC']['nominal']['nugen_e_ME'] = data_path+'L4_IC86.2012_merged_nugen_nue_medium_energy_LE_and_HE_files_49989.npy'
allFiles['MC']['nominal']['nugen_mu_LE'] = data_path+'L4_IC86.2012_merged_nugen_numu_low_energy_LE_and_HE_files_6980.npy'
allFiles['MC']['nominal']['nugen_mu_ME'] = data_path+'L4_IC86.2012_merged_nugen_numu_medium_energy_LE_and_HE_files_49998.npy'
allFiles['MC']['nominal']['nugen_tau_LE'] = data_path+'L4_IC86.2012_merged_nugen_nutau_low_energy_LE_and_HE_files_6969.npy'
allFiles['MC']['nominal']['nugen_tau_ME'] = data_path+'L4_IC86.2012_merged_nugen_nutau_medium_energy_LE_and_HE_files_49997.npy'
allFiles['MC']['nominal']['genie_e'] = data_path+'L4_IC86.2013_merged_genie_LE_and_HE_12550_files_243.npy'
allFiles['MC']['nominal']['genie_mu'] = data_path+'L4_IC86.2013_merged_genie_LE_and_HE_14550_files_614.npy'
allFiles['MC']['nominal']['genie_tau'] = data_path+'L4_IC86.2013_merged_genie_LE_and_HE_16550_files_44.npy'

allFiles['MC']['DomEffUp'] = {}
allFiles['MC']['DomEffUp']['nugen_e_LE'] = data_path+'L4_IC86.2013_merged_nugen_nue_108_low_energy_LE_and_HE_files_6970.npy'
allFiles['MC']['DomEffUp']['nugen_e_ME'] = data_path+'L4_IC86.2013_merged_nugen_nue_108_medium_energy_LE_and_HE_files_30994.npy'
allFiles['MC']['DomEffUp']['nugen_mu_LE'] = data_path+'L4_IC86.2013_merged_nugen_numu_108_low_energy_LE_and_HE_files_6981.npy'
allFiles['MC']['DomEffUp']['nugen_mu_ME'] = data_path+'L4_IC86.2013_merged_nugen_numu_108_medium_energy_LE_and_HE_files_49998.npy'
allFiles['MC']['DomEffUp']['nugen_tau_LE'] = data_path+'L4_IC86.2013_merged_nugen_nutau_108_low_energy_LE_and_HE_files_6977.npy'
allFiles['MC']['DomEffUp']['nugen_tau_ME'] = data_path+'L4_IC86.2013_merged_nugen_nutau_108_medium_energy_LE_and_HE_files_49998.npy'
allFiles['MC']['DomEffUp']['genie_e'] = data_path+'L4_IC86.2013_merged_genie_LE_and_HE_12555_files_217.npy'
allFiles['MC']['DomEffUp']['genie_mu'] = data_path+'L4_IC86.2013_merged_genie_LE_and_HE_14555_files_611.npy'
allFiles['MC']['DomEffUp']['genie_tau'] = data_path+'L4_IC86.2013_merged_genie_LE_and_HE_16555_files_37.npy'

allFiles['MC']['DomEffDown'] = {}
allFiles['MC']['DomEffDown']['nugen_e_LE'] = data_path+'L4_IC86.2013_merged_nugen_nue_090_low_energy_LE_and_HE_files_6952.npy'
allFiles['MC']['DomEffDown']['nugen_e_ME'] = data_path+'L4_IC86.2013_merged_nugen_nue_090_medium_energy_LE_and_HE_files_30994.npy'
allFiles['MC']['DomEffDown']['nugen_mu_LE'] = data_path+'L4_IC86.2013_merged_nugen_numu_090_low_energy_LE_and_HE_files_6983.npy'
allFiles['MC']['DomEffDown']['nugen_mu_ME'] = data_path+'L4_IC86.2013_merged_nugen_numu_090_medium_energy_LE_and_HE_files_30998.npy'
allFiles['MC']['DomEffDown']['nugen_tau_LE'] = data_path+'L4_IC86.2013_merged_nugen_nutau_090_low_energy_LE_and_HE_files_6965.npy'
allFiles['MC']['DomEffDown']['nugen_tau_ME'] = data_path+'L4_IC86.2013_merged_nugen_nutau_090_medium_energy_LE_and_HE_files_30998.npy'
allFiles['MC']['DomEffDown']['genie_e'] = data_path+'L4_IC86.2013_merged_genie_LE_and_HE_12552_files_246.npy'
allFiles['MC']['DomEffDown']['genie_mu'] = data_path+'L4_IC86.2013_merged_genie_LE_and_HE_14552_files_642.npy'
allFiles['MC']['DomEffDown']['genie_tau'] = data_path+'L4_IC86.2013_merged_genie_LE_and_HE_16552_files_48.npy'

allFiles['MC']['Ice_HoleIce100_ScatAbs-7'] = {}
allFiles['MC']['Ice_HoleIce100_ScatAbs-7']['nugen_e_LE'] = data_path+'L4_IC86.2013_merged_nugen_nue_093scatt_abs_low_energy_LE_and_HE_files_6970.npy'
allFiles['MC']['Ice_HoleIce100_ScatAbs-7']['nugen_e_ME'] = data_path+'L4_IC86.2013_merged_nugen_nue_093scatt_abs_medium_energy_LE_and_HE_files_6989.npy'
allFiles['MC']['Ice_HoleIce100_ScatAbs-7']['nugen_mu_LE'] = data_path+'L4_IC86.2013_merged_nugen_numu_093scatt_abs_low_energy_LE_and_HE_files_6980.npy'
allFiles['MC']['Ice_HoleIce100_ScatAbs-7']['nugen_mu_ME'] = data_path+'L4_IC86.2013_merged_nugen_numu_093scatt_abs_medium_energy_LE_and_HE_files_6998.npy'
allFiles['MC']['Ice_HoleIce100_ScatAbs-7']['nugen_tau_LE'] = data_path+'L4_IC86.2013_merged_nugen_nutau_093scatt_abs_low_energy_LE_and_HE_files_6972.npy'
allFiles['MC']['Ice_HoleIce100_ScatAbs-7']['nugen_tau_ME'] = data_path+'L4_IC86.2013_merged_nugen_nutau_093scatt_abs_medium_energy_LE_and_HE_files_6996.npy'
allFiles['MC']['Ice_HoleIce100_ScatAbs-7']['genie_e'] = data_path+'L4_IC86.2013_merged_genie_LE_and_HE_12560_files_252.npy'
allFiles['MC']['Ice_HoleIce100_ScatAbs-7']['genie_mu'] = data_path+'L4_IC86.2013_merged_genie_LE_and_HE_14560_files_618.npy'
allFiles['MC']['Ice_HoleIce100_ScatAbs-7']['genie_tau'] = data_path+'L4_IC86.2013_merged_genie_LE_and_HE_16560_files_41.npy'

allFiles['MC']['Ice_HoleIce100_Scat+10'] = {}
allFiles['MC']['Ice_HoleIce100_Scat+10']['nugen_e_LE'] = data_path+'L4_IC86.2013_merged_nugen_nue_110scatt_low_energy_LE_and_HE_files_6952.npy'
allFiles['MC']['Ice_HoleIce100_Scat+10']['nugen_e_ME'] = data_path+'L4_IC86.2013_merged_nugen_nue_110scatt_medium_energy_LE_and_HE_files_6992.npy'
allFiles['MC']['Ice_HoleIce100_Scat+10']['nugen_mu_LE'] = data_path+ 'L4_IC86.2013_merged_nugen_numu_110scatt_low_energy_LE_and_HE_files_6984.npy'
allFiles['MC']['Ice_HoleIce100_Scat+10']['nugen_mu_ME'] = data_path+'L4_IC86.2013_merged_nugen_numu_110scatt_medium_energy_LE_and_HE_files_6999.npy'
allFiles['MC']['Ice_HoleIce100_Scat+10']['nugen_tau_LE'] = data_path+'L4_IC86.2013_merged_nugen_nutau_110scatt_low_energy_LE_and_HE_files_6968.npy'
allFiles['MC']['Ice_HoleIce100_Scat+10']['nugen_tau_ME'] = data_path+'L4_IC86.2013_merged_nugen_nutau_110scatt_medium_energy_LE_and_HE_files_6996.npy'
allFiles['MC']['Ice_HoleIce100_Scat+10']['genie_e'] = data_path+'L4_IC86.2013_merged_genie_LE_and_HE_12560_files_252.npy'
allFiles['MC']['Ice_HoleIce100_Scat+10']['genie_mu'] = data_path+'L4_IC86.2013_merged_genie_LE_and_HE_14560_files_618.npy'
allFiles['MC']['Ice_HoleIce100_Scat+10']['genie_tau'] = data_path+'L4_IC86.2013_merged_genie_LE_and_HE_16560_files_41.npy'

allFiles['MC']['Ice_HoleIce100_Abs+10'] = {}
allFiles['MC']['Ice_HoleIce100_Abs+10']['nugen_e_LE'] = data_path+'L4_IC86.2013_merged_nugen_nue_110abs_low_energy_LE_and_HE_files_6954.npy'
allFiles['MC']['Ice_HoleIce100_Abs+10']['nugen_e_ME'] = data_path+'L4_IC86.2013_merged_nugen_nue_110abs_medium_energy_LE_and_HE_files_6993.npy'
allFiles['MC']['Ice_HoleIce100_Abs+10']['nugen_mu_LE'] = data_path+'L4_IC86.2013_merged_nugen_numu_110abs_low_energy_LE_and_HE_files_6986.npy'
allFiles['MC']['Ice_HoleIce100_Abs+10']['nugen_mu_ME'] = data_path+'L4_IC86.2013_merged_nugen_numu_110abs_medium_energy_LE_and_HE_files_6999.npy'
allFiles['MC']['Ice_HoleIce100_Abs+10']['nugen_tau_LE'] = data_path+'L4_IC86.2013_merged_nugen_nutau_110abs_low_energy_LE_and_HE_files_6967.npy'
allFiles['MC']['Ice_HoleIce100_Abs+10']['nugen_tau_ME'] = data_path+'L4_IC86.2013_merged_nugen_nutau_110abs_medium_energy_LE_and_HE_files_6998.npy'
allFiles['MC']['Ice_HoleIce100_Abs+10']['genie_e'] = data_path+'L4_IC86.2013_merged_genie_LE_and_HE_12560_files_252.npy'
allFiles['MC']['Ice_HoleIce100_Abs+10']['genie_mu'] = data_path+'L4_IC86.2013_merged_genie_LE_and_HE_14560_files_618.npy'
allFiles['MC']['Ice_HoleIce100_Abs+10']['genie_tau'] = data_path+'L4_IC86.2013_merged_genie_LE_and_HE_16560_files_41.npy'

allFiles['MC']['Ice_HoleIce30_ScatAbs-7'] = {}
allFiles['MC']['Ice_HoleIce30_ScatAbs-7']['nugen_e_LE'] = data_path+'L4_IC86.2013_merged_nugen_nue_093scatt_abs_low_energy_LE_and_HE_files_6970.npy'
allFiles['MC']['Ice_HoleIce30_ScatAbs-7']['nugen_e_ME'] = data_path+'L4_IC86.2013_merged_nugen_nue_093scatt_abs_medium_energy_LE_and_HE_files_6989.npy'
allFiles['MC']['Ice_HoleIce30_ScatAbs-7']['nugen_mu_LE'] = data_path+'L4_IC86.2013_merged_nugen_numu_093scatt_abs_low_energy_LE_and_HE_files_6980.npy'
allFiles['MC']['Ice_HoleIce30_ScatAbs-7']['nugen_mu_ME'] = data_path+'L4_IC86.2013_merged_nugen_numu_093scatt_abs_medium_energy_LE_and_HE_files_6998.npy'
allFiles['MC']['Ice_HoleIce30_ScatAbs-7']['nugen_tau_LE'] = data_path+'L4_IC86.2013_merged_nugen_nutau_093scatt_abs_low_energy_LE_and_HE_files_6972.npy'
allFiles['MC']['Ice_HoleIce30_ScatAbs-7']['nugen_tau_ME'] = data_path+'L4_IC86.2013_merged_nugen_nutau_093scatt_abs_medium_energy_LE_and_HE_files_6996.npy'
allFiles['MC']['Ice_HoleIce30_ScatAbs-7']['genie_e'] = data_path+'L4_IC86.2013_merged_genie_LE_and_HE_12561_files_252.npy'
allFiles['MC']['Ice_HoleIce30_ScatAbs-7']['genie_mu'] = data_path+'L4_IC86.2013_merged_genie_LE_and_HE_14561_files_618.npy'
allFiles['MC']['Ice_HoleIce30_ScatAbs-7']['genie_tau'] = data_path+'L4_IC86.2013_merged_genie_LE_and_HE_16561_files_38.npy'

allFiles['MC']['Ice_HoleIce30_Scat+10'] = {}
allFiles['MC']['Ice_HoleIce30_Scat+10']['nugen_e_LE'] = data_path+'L4_IC86.2013_merged_nugen_nue_110scatt_low_energy_LE_and_HE_files_6952.npy'
allFiles['MC']['Ice_HoleIce30_Scat+10']['nugen_e_ME'] = data_path+'L4_IC86.2013_merged_nugen_nue_110scatt_medium_energy_LE_and_HE_files_6992.npy'
allFiles['MC']['Ice_HoleIce30_Scat+10']['nugen_mu_LE'] = data_path+'L4_IC86.2013_merged_nugen_numu_110scatt_low_energy_LE_and_HE_files_6984.npy'
allFiles['MC']['Ice_HoleIce30_Scat+10']['nugen_mu_ME'] = data_path+'L4_IC86.2013_merged_nugen_numu_110scatt_medium_energy_LE_and_HE_files_6999.npy'
allFiles['MC']['Ice_HoleIce30_Scat+10']['nugen_tau_LE'] = data_path+'L4_IC86.2013_merged_nugen_nutau_110scatt_low_energy_LE_and_HE_files_6968.npy'
allFiles['MC']['Ice_HoleIce30_Scat+10']['nugen_tau_ME'] = data_path+'L4_IC86.2013_merged_nugen_nutau_110scatt_medium_energy_LE_and_HE_files_6996.npy'
allFiles['MC']['Ice_HoleIce30_Scat+10']['genie_e'] = data_path+'L4_IC86.2013_merged_genie_LE_and_HE_12561_files_252.npy'
allFiles['MC']['Ice_HoleIce30_Scat+10']['genie_mu'] = data_path+'L4_IC86.2013_merged_genie_LE_and_HE_14561_files_618.npy'
allFiles['MC']['Ice_HoleIce30_Scat+10']['genie_tau'] = data_path+'L4_IC86.2013_merged_genie_LE_and_HE_16561_files_38.npy'

allFiles['MC']['Ice_HoleIce30_Abs+10'] = {}
allFiles['MC']['Ice_HoleIce30_Abs+10']['nugen_e_LE'] = data_path+'L4_IC86.2013_merged_nugen_nue_110abs_low_energy_LE_and_HE_files_6954.npy'
allFiles['MC']['Ice_HoleIce30_Abs+10']['nugen_e_ME'] = data_path+'L4_IC86.2013_merged_nugen_nue_110abs_medium_energy_LE_and_HE_files_6993.npy'
allFiles['MC']['Ice_HoleIce30_Abs+10']['nugen_mu_LE'] = data_path+'L4_IC86.2013_merged_nugen_numu_110abs_low_energy_LE_and_HE_files_6986.npy'
allFiles['MC']['Ice_HoleIce30_Abs+10']['nugen_mu_ME'] = data_path+'L4_IC86.2013_merged_nugen_numu_110abs_medium_energy_LE_and_HE_files_6999.npy'
allFiles['MC']['Ice_HoleIce30_Abs+10']['nugen_tau_LE'] = data_path+'L4_IC86.2013_merged_nugen_nutau_110abs_low_energy_LE_and_HE_files_6967.npy'
allFiles['MC']['Ice_HoleIce30_Abs+10']['nugen_tau_ME'] = data_path+'L4_IC86.2013_merged_nugen_nutau_110abs_medium_energy_LE_and_HE_files_6998.npy'
allFiles['MC']['Ice_HoleIce30_Abs+10']['genie_e'] = data_path+'L4_IC86.2013_merged_genie_LE_and_HE_12561_files_252.npy'
allFiles['MC']['Ice_HoleIce30_Abs+10']['genie_mu'] = data_path+'L4_IC86.2013_merged_genie_LE_and_HE_14561_files_618.npy'
allFiles['MC']['Ice_HoleIce30_Abs+10']['genie_tau'] = data_path+'L4_IC86.2013_merged_genie_LE_and_HE_16561_files_38.npy'

allFiles['MC']['nominalGammaUp'] = {}
allFiles['MC']['nominalGammaUp']['nugen_e_LE'] = data_path+'L4_IC86.2012_merged_nugen_nue_low_energy_LE_and_HE_files_6953.npy'
allFiles['MC']['nominalGammaUp']['nugen_e_ME'] = data_path+'L4_IC86.2012_merged_nugen_nue_medium_energy_LE_and_HE_files_49989.npy'
allFiles['MC']['nominalGammaUp']['nugen_mu_LE'] = data_path+'L4_IC86.2012_merged_nugen_numu_low_energy_LE_and_HE_files_6980.npy'
allFiles['MC']['nominalGammaUp']['nugen_mu_ME'] = data_path+'L4_IC86.2012_merged_nugen_numu_medium_energy_LE_and_HE_files_49998.npy'
allFiles['MC']['nominalGammaUp']['nugen_tau_LE'] = data_path+'L4_IC86.2012_merged_nugen_nutau_low_energy_LE_and_HE_files_6969.npy'
allFiles['MC']['nominalGammaUp']['nugen_tau_ME'] = data_path+'L4_IC86.2012_merged_nugen_nutau_medium_energy_LE_and_HE_files_49997.npy'
allFiles['MC']['nominalGammaUp']['genie_e'] = data_path+'L4_IC86.2013_merged_genie_LE_and_HE_12550_files_243.npy'
allFiles['MC']['nominalGammaUp']['genie_mu'] = data_path+'L4_IC86.2013_merged_genie_LE_and_HE_14550_files_614.npy'
allFiles['MC']['nominalGammaUp']['genie_tau'] = data_path+'L4_IC86.2013_merged_genie_LE_and_HE_16550_files_44.npy'

allFiles['MC']['nominalGammaDown'] = {}
allFiles['MC']['nominalGammaDown']['nugen_e_LE'] = data_path+'L4_IC86.2012_merged_nugen_nue_low_energy_LE_and_HE_files_6953.npy'
allFiles['MC']['nominalGammaDown']['nugen_e_ME'] = data_path+'L4_IC86.2012_merged_nugen_nue_medium_energy_LE_and_HE_files_49989.npy'
allFiles['MC']['nominalGammaDown']['nugen_mu_LE'] = data_path+'L4_IC86.2012_merged_nugen_numu_low_energy_LE_and_HE_files_6980.npy'
allFiles['MC']['nominalGammaDown']['nugen_mu_ME'] = data_path+'L4_IC86.2012_merged_nugen_numu_medium_energy_LE_and_HE_files_49998.npy'
allFiles['MC']['nominalGammaDown']['nugen_tau_LE'] = data_path+'L4_IC86.2012_merged_nugen_nutau_low_energy_LE_and_HE_files_6969.npy'
allFiles['MC']['nominalGammaDown']['nugen_tau_ME'] = data_path+'L4_IC86.2012_merged_nugen_nutau_medium_energy_LE_and_HE_files_49997.npy'
allFiles['MC']['nominalGammaDown']['genie_e'] = data_path+'L4_IC86.2013_merged_genie_LE_and_HE_12550_files_243.npy'
allFiles['MC']['nominalGammaDown']['genie_mu'] = data_path+'L4_IC86.2013_merged_genie_LE_and_HE_14550_files_614.npy'
allFiles['MC']['nominalGammaDown']['genie_tau'] = data_path+'L4_IC86.2013_merged_genie_LE_and_HE_16550_files_44.npy'


## data

allFiles['Burnsample'] = {}
allFiles['Burnsample']['2012'] = data_path+'L4_IC86.2012_data_8_Runs_merged_LE_or_HE.npy'
allFiles['Burnsample']['2013'] = data_path+'L4_IC86.2013_data_8_Runs_merged_LE_or_HE.npy'
allFiles['Burnsample']['2014'] = data_path+'L4_IC86.2014_data_8_Runs_merged_LE_or_HE.npy'
allFiles['Burnsample']['2015'] = data_path+'L4_IC86.2015_data_8_Runs_merged_LE_or_HE.npy'
allFiles['Burnsample']['2016'] = data_path+'L4_IC86.2016_data_7_Runs_merged_LE_or_HE.npy'

allFiles['Data'] = {}
allFiles['Data']['2012'] = data_path+'L4_IC86.2012_data_merged_LE_or_HE.npy'
allFiles['Data']['2013'] = data_path+'L4_IC86.2013_data_merged_LE_or_HE.npy'
allFiles['Data']['2014'] = data_path+'L4_IC86.2014_data_merged_LE_or_HE.npy'
allFiles['Data']['2015'] = data_path+'L4_IC86.2015_data_merged_LE_or_HE.npy'
allFiles['Data']['2016'] = data_path+'L4_IC86.2016_data_merged_LE_or_HE.npy'


### nfiles normalization for MC datasets

nfiles = {}
nfiles['corsika'] = 1.
nfiles['nominal'] = {'nugen_e_LE':6953. ,'nugen_e_ME': 49989., 'nugen_mu_LE': 6980., 'nugen_mu_ME': 49998., 'nugen_tau_LE':6969.,
                     'nugen_tau_ME': 49997., 'genie_e':243. , 'genie_mu': 614., 'genie_tau':44. }
nfiles['DomEffUp'] = {'nugen_e_LE':6970. ,'nugen_e_ME': 30994., 'nugen_mu_LE': 6981., 'nugen_mu_ME': 49998., 'nugen_tau_LE':6977.,
                      'nugen_tau_ME': 49998., 'genie_e':217. , 'genie_mu': 611., 'genie_tau':37. }
nfiles['DomEffDown'] = {'nugen_e_LE':6952. ,'nugen_e_ME': 30994., 'nugen_mu_LE': 6983., 'nugen_mu_ME': 30998., 'nugen_tau_LE':6965.,
                        'nugen_tau_ME': 30998., 'genie_e':246. , 'genie_mu': 642., 'genie_tau':48. }
nfiles['Ice_HoleIce100_ScatAbs-7'] = {'nugen_e_LE':6970. ,'nugen_e_ME': 6989., 'nugen_mu_LE': 6980., 'nugen_mu_ME': 6998.,
                                      'nugen_tau_LE':6972., 'nugen_tau_ME': 6996., 'genie_e':252. , 'genie_mu': 618., 'genie_tau':41.}
nfiles['Ice_HoleIce100_Scat+10'] = {'nugen_e_LE':6952. ,'nugen_e_ME': 6992., 'nugen_mu_LE': 6984., 'nugen_mu_ME': 6999.,
                                    'nugen_tau_LE':6968., 'nugen_tau_ME': 6996., 'genie_e':252. , 'genie_mu': 618., 'genie_tau':41.}
nfiles['Ice_HoleIce100_Abs+10'] = {'nugen_e_LE':6954. ,'nugen_e_ME': 6993., 'nugen_mu_LE': 6986., 'nugen_mu_ME': 6999.,
                                   'nugen_tau_LE':6967., 'nugen_tau_ME': 6998., 'genie_e':252. , 'genie_mu': 618., 'genie_tau':41. }
nfiles['Ice_HoleIce30_ScatAbs-7'] = {'nugen_e_LE':6970. ,'nugen_e_ME': 6989., 'nugen_mu_LE': 6980., 'nugen_mu_ME': 6998., 
                                     'nugen_tau_LE':6972., 'nugen_tau_ME': 6996., 'genie_e':252. , 'genie_mu': 618., 'genie_tau':38. }
nfiles['Ice_HoleIce30_Scat+10'] = {'nugen_e_LE':6952. ,'nugen_e_ME': 6992., 'nugen_mu_LE': 6984., 'nugen_mu_ME': 6999., 
                                   'nugen_tau_LE':6968.,'nugen_tau_ME': 6996., 'genie_e':252. , 'genie_mu': 618., 'genie_tau':38. }
nfiles['Ice_HoleIce30_Abs+10'] = {'nugen_e_LE':6954. ,'nugen_e_ME': 6993., 'nugen_mu_LE': 6986., 'nugen_mu_ME': 6999., 
                                  'nugen_tau_LE':6967.,'nugen_tau_ME': 6998., 'genie_e':252. , 'genie_mu': 618., 'genie_tau':38. }
nfiles['nominalGammaUp'] = {'nugen_e_LE':6953. ,'nugen_e_ME': 49989., 'nugen_mu_LE': 6980., 'nugen_mu_ME': 49998., 'nugen_tau_LE':6969.,
                     'nugen_tau_ME': 49997., 'genie_e':243. , 'genie_mu': 614., 'genie_tau':44. }
nfiles['nominalGammaDown'] = {'nugen_e_LE':6953. ,'nugen_e_ME': 49989., 'nugen_mu_LE': 6980., 'nugen_mu_ME': 49998., 
                              'nugen_tau_LE':6969.,'nugen_tau_ME': 49997., 'genie_e':243. , 'genie_mu': 614., 'genie_tau':44. }
