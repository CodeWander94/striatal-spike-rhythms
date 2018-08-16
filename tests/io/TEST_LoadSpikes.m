function tests = TEST_LoadSpikes
tests = functiontests(localfunctions);
end

function BasicDataLoadingTest(testCase)

testCase.TestData.S1 = LoadSpikes([]);

verifyEqual(testCase,CheckTS(testCase.TestData.S1),1);
verifyEqual(testCase,length(testCase.TestData.S1.t),92);

cfg = []; cfg.load_questionable_cells = 1;
testCase.TestData.S2 = LoadSpikes(cfg);

verifyEqual(testCase,CheckTS(testCase.TestData.S2),1);
verifyEqual(testCase,length(testCase.TestData.S2.t),128);
end



function setupOnce(testCase)
testCase.TestData.origPath = pwd;
ws = getenv('WORKSPACE');
testFolder1 = fullfile(ws,'testdata','R050-2014-04-02');
cd(testFolder1)
end

function teardownOnce(testCase)
cd(testCase.TestData.origPath);
end
