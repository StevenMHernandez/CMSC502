accuracy_measures = [
    15, 1824.007384, 1907.410399, 1907.410399;
    16, 1925.920728, 1969.558699, 1969.558699;
    17, 1948.714059, 1969.955500, 1969.955500;
    18, 2009.792361, 2185.931969, 2185.931969;
    19, 2089.410146, 2200.467632, 2200.467632;
    21, 2113.540385, 2305.439236, 2305.439236;
    22, 2122.136262, 2307.660320, 2307.660320;
    23, 2139.050036, 2238.366671, 2238.366671;
    24, 2169.661973, 2260.677125, 2260.677125;
    25, 2312.181097, 2369.344695, 2369.344695;
];

figure(1)
hold on

plot(accuracy_measures(:, 1), accuracy_measures(:, 2))
plot(accuracy_measures(:, 1), accuracy_measures(:, 3))
plot(accuracy_measures(:, 1), accuracy_measures(:, 4))

title("Best TSP solution for each Method")
legend("Sequential","Threaded (4 threads)","MPI (4 processes)");

xlabel("# of cities");
ylabel("Distance of best TSP solution");

hold off




figure(2)
hold on

err = ((accuracy_measures(:, 3) - accuracy_measures(:, 2)) ./ accuracy_measures(:, 2)) .* 100;
plot(accuracy_measures(:, 1), err);

title("Percent Error for Threaded Version (4 threads) vs. Sequential Perfect Solution")
xlabel("# of cities");
ylabel("% error");

hold off
