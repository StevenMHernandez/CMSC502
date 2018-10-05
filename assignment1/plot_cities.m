% % sequential
% 
% path = [420.490000 198.710000;319.680000 263.330000;240.010000 315.360000;392.090000 399.720000;476.230000 458.310000;248.060000 486.460000;169.270000 384.690000;148.030000 386.250000;70.770000 402.580000;98.130000 332.460000;72.950000 304.470000;140.690000 278.100000;183.980000 257.920000;143.450000 177.850000;80.450000 201.970000;10.610000 123.340000;67.070000 56.630000;264.180000 45.310000;456.040000 100.780000;445.390000 176.070000;];
% 
% figure(1)
% hold on
% plot(path(:,1), path(:,2), 'o')
% plot(path(:,1), path(:,2))


% threaded
seq_paths = java.util.HashMap;

seq_paths.put(0, [10.610000 123.340000;28.340000 81.010000;49.010000 69.610000;53.830000 65.220000;34.430000 12.460000;67.070000 56.630000;90.160000 121.930000;37.280000 104.320000;]);
seq_paths.put(1, [80.450000 201.970000;57.860000 182.400000;95.800000 139.930000;83.660000 197.370000;96.130000 219.660000;]);
seq_paths.put(2, [72.950000 304.470000;113.770000 326.440000;98.130000 332.460000;84.140000 373.170000;78.310000 366.740000;20.120000 374.530000;39.400000 317.450000;28.630000 261.980000;]);
seq_paths.put(3, [70.770000 402.580000;64.920000 397.250000;4.570000 461.730000;37.200000 474.790000;39.580000 475.180000;121.040000 485.390000;116.410000 446.950000;86.880000 453.630000;75.960000 440.830000;84.070000 415.430000;]);
seq_paths.put(4, [230.210000 33.890000;242.540000 109.870000;224.900000 114.990000;166.070000 117.640000;148.840000 118.050000;166.500000 96.900000;175.190000 94.350000;218.670000 4.110000;]);
seq_paths.put(5, [143.450000 177.850000;242.260000 154.220000;]);
seq_paths.put(6, [140.690000 278.100000;135.170000 271.030000;141.820000 274.190000;183.980000 257.920000;181.150000 277.360000;237.060000 297.290000;240.010000 315.360000;208.460000 348.750000;176.800000 344.120000;143.940000 369.920000;174.110000 300.240000;]);
seq_paths.put(7, [169.270000 384.690000;148.030000 386.250000;169.830000 424.220000;132.100000 439.000000;153.890000 496.130000;248.060000 486.460000;220.220000 466.090000;221.180000 462.170000;201.610000 446.040000;221.450000 440.340000;232.060000 410.290000;248.980000 380.840000;200.720000 407.850000;189.170000 380.720000;]);
seq_paths.put(8, [264.180000 45.310000;267.390000 46.100000;297.960000 92.240000;293.280000 124.100000;315.670000 116.500000;325.730000 125.900000;370.870000 103.100000;356.760000 93.320000;344.710000 85.070000;360.340000 58.860000;344.480000 52.070000;343.850000 49.140000;311.860000 22.830000;266.970000 22.040000;261.300000 41.420000;]);
seq_paths.put(9, [307.290000 149.780000;310.750000 142.330000;327.430000 130.490000;334.110000 146.170000;351.060000 160.140000;366.040000 165.870000;375.510000 185.910000;343.010000 193.140000;321.440000 217.400000;338.550000 242.770000;334.270000 249.890000;290.770000 227.660000;279.330000 209.710000;320.890000 178.640000;]);
seq_paths.put(10, [318.770000 359.360000;337.780000 319.730000;374.540000 315.380000;319.680000 263.330000;289.400000 297.620000;305.530000 314.510000;257.490000 334.690000;287.730000 378.400000;]);
seq_paths.put(11, [257.680000 419.960000;289.540000 439.110000;256.570000 485.800000;329.140000 483.780000;342.900000 455.710000;329.510000 429.690000;320.630000 380.470000;266.580000 379.250000;264.560000 385.530000;]);
seq_paths.put(12, [456.040000 100.780000;499.460000 111.080000;500.000000 104.150000;498.910000 29.390000;492.220000 4.280000;444.750000 40.800000;435.590000 38.480000;445.250000 64.920000;417.040000 118.860000;]);
seq_paths.put(13, [420.490000 198.710000;393.540000 155.460000;388.440000 145.970000;414.130000 167.410000;415.030000 166.840000;445.390000 176.070000;472.300000 226.830000;]);
seq_paths.put(14, [465.580000 361.170000;399.400000 367.000000;460.650000 343.170000;478.340000 295.350000;489.770000 372.550000;]);
seq_paths.put(15, [392.090000 399.720000;404.340000 459.720000;417.040000 462.870000;451.920000 491.840000;492.410000 467.660000;479.420000 459.670000;475.250000 460.260000;476.230000 458.310000;451.350000 425.830000;456.730000 410.300000;436.950000 415.940000;]);



connector = java.util.HashMap;

connector.put(0, [10.610000 123.340000;57.860000 182.400000]);
connector.put(1, [80.450000 201.970000;143.450000 177.850000]);
connector.put(2, [242.260000 154.220000;224.900000 114.990000]);
connector.put(3, [242.540000 109.870000;293.280000 124.100000]);
connector.put(4, [315.670000 116.500000;327.430000 130.490000]);
connector.put(5, [310.750000 142.330000;388.440000 145.970000]);
connector.put(6, [414.130000 167.410000;417.040000 118.860000]);
connector.put(7, [445.250000 64.920000;478.340000 295.350000]);
connector.put(8, [489.770000 372.550000;456.730000 410.300000]);
connector.put(9, [451.350000 425.830000;342.900000 455.710000]);
connector.put(10, [329.510000 429.690000;287.730000 378.400000]);
connector.put(11, [257.490000 334.690000;240.010000 315.360000]);
connector.put(12, [208.460000 348.750000;189.170000 380.720000]);
connector.put(13, [200.720000 407.850000;116.410000 446.950000]);
connector.put(14, [86.880000 453.630000;84.140000 373.170000]);
connector.put(15, [98.130000 332.460000;90.160000 121.930000]);

figure(2)
hold on

for i = [0:15]
    x = seq_paths.get(i);
    x = vertcat(x,x(1,:));
    if ~isempty(x)
        plot(x(:,1), x(:,2), 'o')
        plot(x(:,1), x(:,2))
    end
end
for i = [0:15]
    x = connector.get(i);
    if ~isempty(x)
        plot(x(:,1), x(:,2), 'o')
        plot(x(:,1), x(:,2), '-.')
    end
end

hold off
