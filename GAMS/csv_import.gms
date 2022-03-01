$call csv2gdx export_tax.csv      id=param_exports_Tax_Table    index=1     value=2..lastCol useHeader=y output=gdx/param_exports_Tax_Table.gdx
$call csv2gdx hh_endowments.csv   id=param_endowment_Table      index=1     value=2..lastCol useHeader=y output=gdx/param_endowment_Table.gdx
$call csv2gdx hh_tax_rates.csv    id=param_hh_tax_rates_Table   index=1     value=2..lastCol useHeader=y output=gdx/param_hh_tax_rates_Table.gdx
$call csv2gdx IO_Export.csv       id=param_exports_Table        index=1     value=2..lastCol useHeader=y output=gdx/param_exports_Table.gdx
$call csv2gdx IO_FactorUse.csv    id=param_factor_use_Table     index=1     value=2..lastCol useHeader=y output=gdx/param_factor_use_Table.gdx
$call csv2gdx IO_hh_pfce.csv      id=param_pfce_Table           index=1,2   value=3..lastCol useHeader=y output=gdx/param_pfce_Table.gdx valueDim=y
$call csv2gdx IO_TaxRate.csv      id=param_Tax_Rate_Table       index=1,2   value=3..lastCol useHeader=y output=gdx/param_Tax_Rate_Table.gdx  valueDim=y
$call csv2gdx IO_IUseData.csv     id=param_IO_Table             index=1,2   value=3..lastCol useHeader=y output=gdx/param_IO_Table.gdx valueDim=y
$call csv2gdx sigma_dm_fu.csv     id=param_sigma_dm_fu_Table    index=1     value=2..lastCol useHeader=y output=gdx/param_sigma_dm_fu_Table.gdx
$call csv2gdx sigma_dm_io.csv     id=param_sigma_dm_io_Table    index=1     value=2..lastCol useHeader=y output=gdx/param_sigma_dm_io_Table.gdx
$call csv2gdx sigma_iu_va.csv     id=param_sigma_iu_va_Table    index=1     value=2..lastCol useHeader=y output=gdx/param_sigma_iu_va_Table.gdx
$call csv2gdx sigma_kl.csv        id=param_sigma_kl_Table       index=1     value=2..lastCol useHeader=y output=gdx/param_sigma_kl_Table.gdx