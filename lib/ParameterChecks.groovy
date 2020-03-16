class ParameterChecks {
	static void checkParams(params) {
		assert params.input, "Please specify path to input with --input parameter"
        assert params.samples, "Please specify path to samples.tsv with --samples parameter"
        assert !(params.Gmodel && !params.SNPs), "Please specify path to variants with --SNPs parameter when running --Gmodel"
        assert !(params.GxE && !params.SNPs), "Please specify path to variants with --SNPs parameter when running --GxE"
        assert !(params.noCpG && params.noCHG && params.noCHH), "Please specify at least one methylation context for analysis!"
        assert params.kplots instanceof Integer && params.kplots >= 0, "Number of kplots must be a non-negative integer!"
        assert params.coverage instanceof Integer && params.coverage >= 0, "Input coverage filter must be a non-negative integer!"
        assert params.distance instanceof Integer && params.distance >= 0, "Distance parameter must be a non-negative integer!"
        assert Double.valueOf(params.proportion) >= 0.0d && Double.valueOf(params.proportion) <= 1.0d, "Proportion of shared regions must be a decimal in the range of 0 and 1!"
        assert Double.valueOf(params.input_FDR) >= 0.0d && Double.valueOf(params.input_FDR) <= 1.0d, "Input FDR filter must be a decimal in the range of 0 and 1!"
        assert Double.valueOf(params.output_FDR) >= 0.0d && Double.valueOf(params.output_FDR) <= 1.0d, "Output FDR filter must be a decimal in the range of 0 and 1!"
        assert Double.valueOf(params.Emodel_pv) >= 0.0d && Double.valueOf(params.Emodel_pv) <= 1.0d, "Emodel p-value must be a decimal in the range of 0 and 1!"
        assert Double.valueOf(params.Gmodel_pv) >= 0.0d && Double.valueOf(params.Gmodel_pv) <= 1.0d, "Gmodel p-value must be a decimal in the range of 0 and 1!"
        assert Double.valueOf(params.GxE_pv) >= 0.0d && Double.valueOf(params.GxE_pv) <= 1.0d, "GxE p-value must be a decimal in the range of 0 and 1!"
        assert Double.valueOf(params.max_missing) >= 0.0d && Double.valueOf(params.max_missing) <= 1.0d, "--max_missing parameter must be a decimal in the range of 0 and 1!"
        assert params.mac instanceof Integer && params.mac >= 0, "--mac parameter must be a non-negative integer!"
        assert params.minQ instanceof Integer && params.minQ >= 0, "--minQ parameter must be a non-negative integer!"
        assert params.take instanceof Integer && params.take >= 0, "--take parameter must be a non-negative integer!"
        assert params.fork instanceof Integer && params.fork >= 0, "--fork parameter must be a non-negative integer!"
	}
}