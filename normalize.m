function wave = normalize(wave,N0,dx)

wave = wave.*sqrt(N0./(dx.*sum(abs(wave).^2)));
