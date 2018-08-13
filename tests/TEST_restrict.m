% setup
empty_ts = ts;

one_ts_data = (0:1:10)';
one_ts = ts; one_ts.t{1} = one_ts_data; one_ts.label{1} = 'label1';

rvec = [0 10]; rvec_iv = iv(rvec);

%% test if format is preserved with raw input times
empty_ts_out = restrict(empty_ts,rvec(1),rvec(2));
assert(CheckTS(empty_ts_out) == 1);

%% test if format is preserved with iv input times
empty_ts_out = restrict(empty_ts,rvec_iv);
assert(CheckTS(empty_ts_out) == 1);
