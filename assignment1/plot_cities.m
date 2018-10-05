% sequential

path = [420.490000 198.710000;319.680000 263.330000;240.010000 315.360000;392.090000 399.720000;476.230000 458.310000;248.060000 486.460000;169.270000 384.690000;148.030000 386.250000;70.770000 402.580000;98.130000 332.460000;72.950000 304.470000;140.690000 278.100000;183.980000 257.920000;143.450000 177.850000;80.450000 201.970000;10.610000 123.340000;67.070000 56.630000;264.180000 45.310000;456.040000 100.780000;445.390000 176.070000;];

figure(1)
hold on
plot(path(:,1), path(:,2), 'o')
plot(path(:,1), path(:,2))


% threaded
seq_paths = java.util.HashMap;

seq_paths.put(0, [10.610000 123.340000;37.280000 104.320000;90.160000 121.930000;91.160000 112.930000;67.070000 56.630000;77.740000 38.590000;78.560000 18.860000;34.430000 12.460000;3.060000 5.340000;28.130000 29.080000;37.370000 51.670000;53.830000 65.220000;49.010000 69.610000;28.340000 81.010000;]);
seq_paths.put(1, [80.450000 201.970000;96.130000 219.660000;107.410000 238.260000;80.680000 252.970000;26.960000 225.930000;52.030000 192.990000;57.860000 182.400000;71.270000 181.820000;92.630000 172.320000;83.470000 138.770000;95.800000 139.930000;124.050000 164.690000;117.420000 194.410000;102.180000 202.090000;83.660000 197.370000;]);
seq_paths.put(2, [72.950000 304.470000;61.760000 295.840000;28.630000 261.980000;39.400000 317.450000;56.150000 324.490000;57.090000 338.120000;62.770000 343.680000;34.120000 343.650000;20.120000 374.530000;78.310000 366.740000;84.140000 373.170000;92.260000 342.880000;98.130000 332.460000;113.770000 326.440000;101.020000 307.190000;95.840000 310.430000;]);
seq_paths.put(3, [70.770000 402.580000;64.920000 397.250000;4.570000 461.730000;37.200000 474.790000;39.580000 475.180000;75.960000 440.830000;86.880000 453.630000;121.040000 485.390000;116.410000 446.950000;118.740000 415.330000;84.070000 415.430000;]);
seq_paths.put(4, [230.210000 33.890000;218.670000 4.110000;206.870000 5.190000;167.360000 47.850000;148.600000 59.550000;175.190000 94.350000;166.500000 96.900000;155.810000 112.450000;148.840000 118.050000;166.070000 117.640000;224.900000 114.990000;242.160000 114.680000;242.540000 109.870000;239.600000 86.990000;244.990000 80.740000;250.520000 62.830000;]);
seq_paths.put(5, [143.450000 177.850000;154.360000 132.630000;183.390000 145.910000;242.260000 154.220000;249.250000 147.300000;249.300000 169.150000;191.660000 192.630000;159.540000 218.420000;]);
seq_paths.put(6, [140.690000 278.100000;132.070000 325.710000;152.630000 329.420000;143.940000 369.920000;154.440000 369.360000;182.670000 360.140000;176.800000 344.120000;198.680000 330.270000;208.460000 348.750000;240.010000 315.360000;237.060000 297.290000;250.330000 289.160000;183.980000 257.920000;181.150000 277.360000;174.110000 300.240000;141.820000 274.190000;135.170000 271.030000;]);
seq_paths.put(7, [169.270000 384.690000;148.030000 386.250000;130.590000 389.680000;173.760000 407.920000;169.830000 424.220000;153.880000 444.640000;132.100000 439.000000;150.400000 452.700000;153.890000 496.130000;164.090000 466.120000;180.200000 439.440000;201.610000 446.040000;215.100000 467.410000;220.220000 466.090000;248.060000 486.460000;221.180000 462.170000;221.450000 440.340000;233.170000 424.850000;232.060000 410.290000;248.980000 380.840000;230.210000 390.980000;200.720000 407.850000;193.460000 387.700000;189.170000 380.720000;]);
seq_paths.put(8, [264.180000 45.310000;267.390000 46.100000;267.590000 79.480000;260.060000 105.900000;297.960000 92.240000;293.280000 124.100000;315.670000 116.500000;325.730000 125.900000;370.870000 103.100000;356.760000 93.320000;344.710000 85.070000;364.460000 71.680000;360.340000 58.860000;344.480000 52.070000;343.850000 49.140000;330.430000 20.570000;323.830000 5.490000;311.860000 22.830000;278.810000 9.750000;266.970000 22.040000;261.300000 41.420000;]);
seq_paths.put(9, [307.290000 149.780000;310.750000 142.330000;292.830000 134.570000;305.420000 131.310000;327.430000 130.490000;334.110000 146.170000;349.680000 157.700000;351.060000 160.140000;366.040000 165.870000;375.510000 185.910000;343.010000 193.140000;321.440000 217.400000;346.810000 222.670000;338.550000 242.770000;334.270000 249.890000;302.550000 247.480000;290.770000 227.660000;279.890000 214.530000;279.330000 209.710000;260.460000 201.850000;298.540000 182.550000;320.890000 178.640000;327.520000 158.330000;]);
seq_paths.put(10, [318.770000 359.360000;257.490000 334.690000;310.180000 322.240000;305.530000 314.510000;310.420000 299.110000;289.400000 297.620000;290.370000 266.120000;305.840000 264.060000;319.680000 263.330000;350.900000 285.000000;374.540000 315.380000;365.360000 320.230000;337.780000 319.730000;]);
seq_paths.put(11, [257.680000 419.960000;264.560000 385.530000;266.580000 379.250000;287.730000 378.400000;320.630000 380.470000;330.230000 381.480000;343.780000 405.370000;313.920000 418.460000;328.550000 429.130000;329.510000 429.690000;334.500000 439.810000;342.900000 455.710000;329.140000 483.780000;256.570000 485.800000;277.280000 460.000000;289.540000 439.110000;267.650000 422.350000;]);
seq_paths.put(12, [456.040000 100.780000;462.050000 86.900000;499.460000 111.080000;500.000000 104.150000;498.910000 29.390000;492.220000 4.280000;444.750000 40.800000;435.590000 38.480000;432.930000 49.220000;445.250000 64.920000;449.040000 70.200000;405.020000 68.020000;417.040000 118.860000;431.300000 103.770000;442.790000 95.170000;]);
seq_paths.put(13, [420.490000 198.710000;415.390000 198.710000;414.130000 167.410000;415.030000 166.840000;393.540000 155.460000;388.440000 145.970000;413.820000 147.190000;445.390000 176.070000;492.460000 170.780000;472.300000 226.830000;432.630000 247.470000;]);
seq_paths.put(14, [465.580000 361.170000;460.650000 343.170000;433.100000 361.200000;419.470000 362.820000;399.400000 367.000000;391.680000 360.430000;383.020000 350.290000;401.440000 339.710000;397.340000 275.150000;441.810000 313.770000;454.470000 311.990000;478.340000 295.350000;493.270000 354.680000;489.770000 372.550000;]);
seq_paths.put(15, [392.090000 399.720000;404.340000 459.720000;417.040000 462.870000;451.920000 491.840000;492.410000 467.660000;479.420000 459.670000;476.230000 458.310000;475.250000 460.260000;471.670000 458.350000;455.050000 437.300000;451.350000 425.830000;456.730000 410.300000;436.950000 415.940000;418.820000 409.520000;392.180000 387.760000;]);

figure(2)

for i = [0:15]
    x = seq_paths.get(i);
    if ~isempty(x)
        hold on
        plot(x(:,1), x(:,2), 'o')
        hold on
        plot(x(:,1), x(:,2))
    end
end

hold off