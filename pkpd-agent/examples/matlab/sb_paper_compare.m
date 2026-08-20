function sb_paper_compare(projPath, vpopXlsx, outDir)
%SB_PAPER_COMPARE Reproduce the paper's Fig 4 cascade at FULL Vpop and dump the
%   per-patient response of each arm, so my pipeline's model output can be diffed
%   number-for-number against the published model-prediction bars (Fig 5 / Fig 6).
%
%   The paper (Bedathuru et al., npj Syst Biol Appl 2024) simulates one Vpop
%   sequentially (their Fig 4):
%       All comers (n=300)
%         --MTX-->  MTX-IR  = ACR<50 at Wk12                    (Fig 5A: MTX naive)
%         --ADA-->  ADA-IR  = DAS28-CRP>3.2 & ACR<50 at Wk24    (Fig 5B: ADA naive)
%       TCZ on MTX-IR                                            (Fig 5C, n=251)
%       TCZ on (MTX-IR & ADA-IR)                                (Fig 6 validation, n=216)
%
%   This script runs three FULL-population passes (no subsampling) and writes three
%   per-patient CSVs. All population set-logic (the IR masks and their intersection)
%   is done afterwards by paper_compare.py, which needs no MATLAB. The three arms:
%       1. MTX_15mg_Q1W_SC_t200            readout Wk12 (day 284)
%       2. ADA40mg_Q2W_SC_t200             readout Wk24 (day 368)
%       3. MTX + TCZ@day285 (sequential)   readout Wk24-post-switch (day 453)
%   Arm 3 doses the WHOLE Vpop with MTX then TCZ; the ACR20/50/70 columns are the
%   model's CONTINUOUS response flags read at day 453, i.e. the post-TCZ response of
%   EVERY patient. The Python step then selects the paper's IR sub-populations from
%   arms 1-2 and reads arm 3's response over them. Nothing here is placebo-corrected:
%   these are raw model outputs, to be compared to the BLACK (Simulation) bars.
%
%   Usage (Windows, pkpd-ml conda env, from the examples\matlab folder):
%       sb_paper_compare('..\..\..\RA-QSP-Model\Vantage RA QSP Model v1.0.sbproj', ...
%                        '..\..\..\RA-QSP-Model\Vpop1.xlsx', '.')
%   then, anywhere:  python -m examples.paper_compare --dir <outDir>

    if nargin < 3 || isempty(outDir), outDir = '.'; end
    if ~exist(outDir, 'dir'), mkdir(outDir); end

    sb_load(projPath);                 % stashes the model in base as 'sbmodel'

    mtxCsv = fullfile(outDir, 'arm_mtx.csv');
    adaCsv = fullfile(outDir, 'arm_ada.csv');
    tczCsv = fullfile(outDir, 'arm_tcz.csv');

    fprintf('\n== ARM 1/3: MTX naive, readout Wk12 (day 284), full Vpop ==\n');
    sb_run_vpop(vpopXlsx, 'MTX_15mg_Q1W_SC_t200', 320, 200, 284, mtxCsv, 0);

    fprintf('\n== ARM 2/3: ADA naive, readout Wk24 (day 368), full Vpop ==\n');
    sb_run_vpop(vpopXlsx, 'ADA40mg_Q2W_SC_t200', 400, 200, 368, adaCsv, 0);

    fprintf('\n== ARM 3/3: MTX then TCZ@day285, readout day 453 (Wk24 post-switch) ==\n');
    sb_run_vpop(vpopXlsx, 'MTX_15mg_Q1W_SC_t200;TCZ8mgkg_Q4W_IV_t200@285', ...
                500, 200, 453, tczCsv, 0);

    fprintf('\n== done ==\n');
    fprintf('wrote:\n  %s\n  %s\n  %s\n', mtxCsv, adaCsv, tczCsv);
    fprintf('next:  python -m examples.paper_compare --dir "%s"\n', outDir);
end
