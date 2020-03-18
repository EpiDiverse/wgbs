class ParameterChecks {
	static void checkParamsLocal(params) {
		assert params.input, "Please specify path to input reads with --input parameter"
        assert params.reference, "Please specify path to reference genome with --reference parameter"
	}
    static void checkParamsEpi(params) {
		assert params.input, "Please specify path to input reads with --input parameter"
        assert params.reference || params.thlaspi || params.fragaria || params.populus, "Please specify path to reference genome with --reference parameter"
    }
    static void checkParams(params) {
        assert !(params.noCpG && params.noCHG && params.noCHH), "Please specify at least one methylation context for analysis!"
        assert params.forward instanceof String, "Please specify a valid forward adapter sequence for trimming"
        assert params.reverse instanceof String, "Please specify a valid reverse adapter sequence for trimming"
        assert params.clip5 instanceof Integer && params.clip5 >= 0, "--clip5 parameter must be a non-negative integer!"
        assert params.clip3 instanceof Integer && params.clip3 >= 0, "--clip3 parameter must be a non-negative integer!"
        assert params.minOver instanceof Integer && params.minOver > 0, "--minOver parameter must be a positive integer!"
        assert params.minLeng instanceof Integer && params.minLeng > 0, "--minLeng parameter must be a positive integer!"
        assert params.minQual instanceof Integer && params.minQual >= 0, "--minQual parameter must be a non-negative integer!"
        assert params.minIns instanceof Integer && params.minIns >= 0, "--minIns parameter must be a non-negative integer!"
        assert params.maxIns instanceof Integer && params.maxIns >= 0, "--maxIns parameter must be a non-negative integer!"
        assert params.maxErrors instanceof Integer && (params.maxErrors == -1 || params.maxErrors >= 0), "--maxErrors parameter must be a non-negative integer (or can be disabled with -1)!"
        assert params.minAccuracy instanceof Integer && params.minAccuracy >= 0, "--minAccuracy parameter must be a non-negative integer!"
        assert Double.valueOf(params.XF) >= 0.0d && Double.valueOf(params.XF) <= 1.0d, "--XF proportion must be a decimal in the range of 0 and 1!"
        assert params.take instanceof Integer && params.take > 0, "--take parameter must be a positive integer!"
        assert !params.fork || (params.fork instanceof Integer && params.fork > 0), "--fork parameter must be a positive integer!"
	}
}