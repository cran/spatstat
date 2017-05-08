
/* 
   Native symbol registration table for spatstat package

   Automatically generated - do not edit this file!

*/

#include "proto.h"
#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/*  
   See proto.h for declarations for the native routines registered below.
*/

static const R_CMethodDef CEntries[] = {
    {"acrdenspt",        (DL_FUNC) &acrdenspt,        10},
    {"acrsmoopt",        (DL_FUNC) &acrsmoopt,        10},
    {"adenspt",          (DL_FUNC) &adenspt,           7},
    {"areaBdif",         (DL_FUNC) &areaBdif,         11},
    {"areadifs",         (DL_FUNC) &areadifs,          7},
    {"asmoopt",          (DL_FUNC) &asmoopt,           8},
    {"auctionbf",        (DL_FUNC) &auctionbf,         7},
    {"awtcrdenspt",      (DL_FUNC) &awtcrdenspt,      11},
    {"awtcrsmoopt",      (DL_FUNC) &awtcrsmoopt,      11},
    {"awtdenspt",        (DL_FUNC) &awtdenspt,         8},
    {"awtsmoopt",        (DL_FUNC) &awtsmoopt,         9},
    {"bdrymask",         (DL_FUNC) &bdrymask,          4},
    {"Cbiform",          (DL_FUNC) &Cbiform,           6},
    {"Cclosepaircounts", (DL_FUNC) &Cclosepaircounts,  5},
    {"Ccountends",       (DL_FUNC) &Ccountends,       14},
    {"Ccrossdist",       (DL_FUNC) &Ccrossdist,        8},
    {"Ccrosspaircounts", (DL_FUNC) &Ccrosspaircounts,  8},
    {"CcrossPdist",      (DL_FUNC) &CcrossPdist,      10},
    {"Cidw",             (DL_FUNC) &Cidw,             14},
    {"ClineMquad",       (DL_FUNC) &ClineMquad,       23},
    {"Clinequad",        (DL_FUNC) &Clinequad,        18},
    {"ClineRMquad",      (DL_FUNC) &ClineRMquad,      23},
    {"ClineRquad",       (DL_FUNC) &ClineRquad,       18},
    {"Clixellate",       (DL_FUNC) &Clixellate,       16},
    {"cocoGraph",        (DL_FUNC) &cocoGraph,         6},
    {"cocoImage",        (DL_FUNC) &cocoImage,         3},
    {"Corput",           (DL_FUNC) &Corput,            3},
    {"Cpairdist",        (DL_FUNC) &Cpairdist,         5},
    {"CpairPdist",       (DL_FUNC) &CpairPdist,        7},
    {"Cquadform",        (DL_FUNC) &Cquadform,         5},
    {"crdenspt",         (DL_FUNC) &crdenspt,          9},
    {"crosscount",       (DL_FUNC) &crosscount,        8},
    {"crsmoopt",         (DL_FUNC) &crsmoopt,         10},
    {"CspaSumSymOut",    (DL_FUNC) &CspaSumSymOut,     9},
    {"CspaWtSumSymOut",  (DL_FUNC) &CspaWtSumSymOut,  13},
    {"Csum2outer",       (DL_FUNC) &Csum2outer,        6},
    {"Csumouter",        (DL_FUNC) &Csumouter,         4},
    {"Csumsymouter",     (DL_FUNC) &Csumsymouter,      4},
    {"Cwsum2outer",      (DL_FUNC) &Cwsum2outer,       7},
    {"Cwsumouter",       (DL_FUNC) &Cwsumouter,        5},
    {"Cwsumsymouter",    (DL_FUNC) &Cwsumsymouter,     5},
    {"Cxypolyselfint",   (DL_FUNC) &Cxypolyselfint,   11},
    {"D3crossdist",      (DL_FUNC) &D3crossdist,      10},
    {"D3crossPdist",     (DL_FUNC) &D3crossPdist,     13},
    {"D3pairdist",       (DL_FUNC) &D3pairdist,        6},
    {"D3pairPdist",      (DL_FUNC) &D3pairPdist,       9},
    {"Ddist2dpath",      (DL_FUNC) &Ddist2dpath,       7},
    {"delta2area",       (DL_FUNC) &delta2area,       10},
    {"denspt",           (DL_FUNC) &denspt,            6},
    {"digberJ",          (DL_FUNC) &digberJ,           6},
    {"dinfty_R",         (DL_FUNC) &dinfty_R,          3},
    {"discareapoly",     (DL_FUNC) &discareapoly,     12},
    {"discs2grid",       (DL_FUNC) &discs2grid,       11},
    {"distmapbin",       (DL_FUNC) &distmapbin,        9},
    {"dwpure",           (DL_FUNC) &dwpure,            6},
    {"Ediggatsti",       (DL_FUNC) &Ediggatsti,       10},
    {"Ediggra",          (DL_FUNC) &Ediggra,          11},
    {"Efiksel",          (DL_FUNC) &Efiksel,           9},
    {"Egeyer",           (DL_FUNC) &Egeyer,           11},
    {"exact_dt_R",       (DL_FUNC) &exact_dt_R,       14},
    {"fardist2grid",     (DL_FUNC) &fardist2grid,     10},
    {"fardistgrid",      (DL_FUNC) &fardistgrid,      10},
    {"Fclosepairs",      (DL_FUNC) &Fclosepairs,      16},
    {"Fcrosspairs",      (DL_FUNC) &Fcrosspairs,      19},
    {"Gdenspt",          (DL_FUNC) &Gdenspt,           5},
    {"Gsmoopt",          (DL_FUNC) &Gsmoopt,           7},
    {"Gwtdenspt",        (DL_FUNC) &Gwtdenspt,         6},
    {"Gwtsmoopt",        (DL_FUNC) &Gwtsmoopt,         8},
    {"hasX3close",       (DL_FUNC) &hasX3close,        6},
    {"hasX3pclose",      (DL_FUNC) &hasX3pclose,       7},
    {"hasXclose",        (DL_FUNC) &hasXclose,         5},
    {"hasXpclose",       (DL_FUNC) &hasXpclose,        6},
    {"hasXY3close",      (DL_FUNC) &hasXY3close,      10},
    {"hasXY3pclose",     (DL_FUNC) &hasXY3pclose,     11},
    {"hasXYclose",       (DL_FUNC) &hasXYclose,        8},
    {"hasXYpclose",      (DL_FUNC) &hasXYpclose,       9},
    {"Idist2dpath",      (DL_FUNC) &Idist2dpath,       7},
    {"idwloo",           (DL_FUNC) &idwloo,            8},
    {"KborderD",         (DL_FUNC) &KborderD,          8},
    {"KborderI",         (DL_FUNC) &KborderI,          8},
    {"knnd3D",           (DL_FUNC) &knnd3D,            8},
    {"knndMD",           (DL_FUNC) &knndMD,            6},
    {"knndsort",         (DL_FUNC) &knndsort,          6},
    {"knnGinterface",    (DL_FUNC) &knnGinterface,    15},
    {"knnsort",          (DL_FUNC) &knnsort,           7},
    {"knnw3D",           (DL_FUNC) &knnw3D,            8},
    {"knnwMD",           (DL_FUNC) &knnwMD,            7},
    {"knnX3Dinterface",  (DL_FUNC) &knnX3Dinterface,  17},
    {"knnXinterface",    (DL_FUNC) &knnXinterface,    15},
    {"KnoneD",           (DL_FUNC) &KnoneD,            6},
    {"KnoneI",           (DL_FUNC) &KnoneI,            6},
    {"knownCif",         (DL_FUNC) &knownCif,          2},
    {"KrectDbl",         (DL_FUNC) &KrectDbl,         17},
    {"KrectInt",         (DL_FUNC) &KrectInt,         17},
    {"KrectWtd",         (DL_FUNC) &KrectWtd,         18},
    {"Kwborder",         (DL_FUNC) &Kwborder,          9},
    {"Kwnone",           (DL_FUNC) &Kwnone,            7},
    {"lincrossdist",     (DL_FUNC) &lincrossdist,     16},
    {"linearradius",     (DL_FUNC) &linearradius,      8},
    {"linknncross",      (DL_FUNC) &linknncross,      16},
    {"linknnd",          (DL_FUNC) &linknnd,          13},
    {"linndcross",       (DL_FUNC) &linndcross,       18},
    {"linndxcross",      (DL_FUNC) &linndxcross,      20},
    {"linnndist",        (DL_FUNC) &linnndist,        13},
    {"linnnwhich",       (DL_FUNC) &linnnwhich,       14},
    {"linpairdist",      (DL_FUNC) &linpairdist,      12},
    {"linSnndwhich",     (DL_FUNC) &linSnndwhich,     15},
    {"locpcfx",          (DL_FUNC) &locpcfx,          12},
    {"locprod",          (DL_FUNC) &locprod,           7},
    {"locWpcfx",         (DL_FUNC) &locWpcfx,         13},
    {"locxprod",         (DL_FUNC) &locxprod,         10},
    {"maxnnd2",          (DL_FUNC) &maxnnd2,           5},
    {"maxPnnd2",         (DL_FUNC) &maxPnnd2,          5},
    {"minnnd2",          (DL_FUNC) &minnnd2,           5},
    {"minPnnd2",         (DL_FUNC) &minPnnd2,          5},
    {"nnd3D",            (DL_FUNC) &nnd3D,             7},
    {"nndistsort",       (DL_FUNC) &nndistsort,        5},
    {"nndMD",            (DL_FUNC) &nndMD,             5},
    {"nnGinterface",     (DL_FUNC) &nnGinterface,     14},
    {"nnw3D",            (DL_FUNC) &nnw3D,             7},
    {"nnwhichsort",      (DL_FUNC) &nnwhichsort,       5},
    {"nnwMD",            (DL_FUNC) &nnwMD,             6},
    {"nnX3Dinterface",   (DL_FUNC) &nnX3Dinterface,   16},
    {"nnXinterface",     (DL_FUNC) &nnXinterface,     14},
    {"paircount",        (DL_FUNC) &paircount,         5},
    {"poly2imA",         (DL_FUNC) &poly2imA,          7},
    {"poly2imI",         (DL_FUNC) &poly2imI,          6},
    {"ps_exact_dt_R",    (DL_FUNC) &ps_exact_dt_R,    13},
    {"RcallF3",          (DL_FUNC) &RcallF3,          17},
    {"RcallF3cen",       (DL_FUNC) &RcallF3cen,       20},
    {"RcallG3",          (DL_FUNC) &RcallG3,          17},
    {"RcallG3cen",       (DL_FUNC) &RcallG3cen,       19},
    {"RcallK3",          (DL_FUNC) &RcallK3,          17},
    {"Rcallpcf3",        (DL_FUNC) &Rcallpcf3,        18},
    {"ripleybox",        (DL_FUNC) &ripleybox,        11},
    {"ripleypoly",       (DL_FUNC) &ripleypoly,       11},
    {"scantrans",        (DL_FUNC) &scantrans,        11},
    {"seg2pixI",         (DL_FUNC) &seg2pixI,          8},
    {"seg2pixL",         (DL_FUNC) &seg2pixL,         11},
    {"seg2pixN",         (DL_FUNC) &seg2pixN,          9},
    {"segdens",          (DL_FUNC) &segdens,          10},
    {"smoopt",           (DL_FUNC) &smoopt,            8},
    {"trigraf",          (DL_FUNC) &trigraf,          10},
    {"trigrafS",         (DL_FUNC) &trigrafS,         10},
    {"wtcrdenspt",       (DL_FUNC) &wtcrdenspt,       10},
    {"wtcrsmoopt",       (DL_FUNC) &wtcrsmoopt,       11},
    {"wtdenspt",         (DL_FUNC) &wtdenspt,          7},
    {"wtsmoopt",         (DL_FUNC) &wtsmoopt,          9},
    {"xypsi",            (DL_FUNC) &xypsi,            10},
    {"xysegint",         (DL_FUNC) &xysegint,         16},
    {"xysegXint",        (DL_FUNC) &xysegXint,        11},
    {"xysi",             (DL_FUNC) &xysi,             12},
    {"xysiANY",          (DL_FUNC) &xysiANY,          12},
    {"xysxi",            (DL_FUNC) &xysxi,             7},
    {"Clinvwhichdist",   (DL_FUNC) &Clinvwhichdist,   12},
    {"Clinvdist",        (DL_FUNC) &Clinvdist,        11},
    {NULL, NULL, 0}
};

static const R_CallMethodDef CallEntries[] = {
    {"close3IJpairs",        (DL_FUNC) &close3IJpairs,         5},
    {"close3pairs",          (DL_FUNC) &close3pairs,           5},
    {"cross3IJpairs",        (DL_FUNC) &cross3IJpairs,         8},
    {"cross3pairs",          (DL_FUNC) &cross3pairs,           8},
    {"Cwhist",               (DL_FUNC) &Cwhist,                3},
    {"Cxysegint",            (DL_FUNC) &Cxysegint,             9},
    {"CxysegXint",           (DL_FUNC) &CxysegXint,            5},
    {"graphVees",            (DL_FUNC) &graphVees,             3},
    {"PerfectDGS",           (DL_FUNC) &PerfectDGS,            4},
    {"PerfectDiggleGratton", (DL_FUNC) &PerfectDiggleGratton,  6},
    {"PerfectHardcore",      (DL_FUNC) &PerfectHardcore,       4},
    {"PerfectPenttinen",     (DL_FUNC) &PerfectPenttinen,      5},
    {"PerfectStrauss",       (DL_FUNC) &PerfectStrauss,        5},
    {"PerfectStraussHard",   (DL_FUNC) &PerfectStraussHard,    6},
    {"thinjumpequal",        (DL_FUNC) &thinjumpequal,         3},
    {"triDgraph",            (DL_FUNC) &triDgraph,             4},
    {"triDRgraph",           (DL_FUNC) &triDRgraph,            5},
    {"triograph",            (DL_FUNC) &triograph,             3},
    {"trioxgraph",           (DL_FUNC) &trioxgraph,            4},
    {"VcloseIJDpairs",       (DL_FUNC) &VcloseIJDpairs,        4},
    {"VcloseIJpairs",        (DL_FUNC) &VcloseIJpairs,         4},
    {"Vclosepairs",          (DL_FUNC) &Vclosepairs,           4},
    {"Vclosethresh",         (DL_FUNC) &Vclosethresh,          5},
    {"VcrossIJDpairs",       (DL_FUNC) &VcrossIJDpairs,        6},
    {"VcrossIJpairs",        (DL_FUNC) &VcrossIJpairs,         6},
    {"Vcrosspairs",          (DL_FUNC) &Vcrosspairs,           6},
    {"xmethas",              (DL_FUNC) &xmethas,              25},
    {NULL, NULL, 0}
};

void R_init_spatstat(DllInfo *dll)
{
    R_registerRoutines(dll, CEntries, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
