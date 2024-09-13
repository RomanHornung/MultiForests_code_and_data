// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include "../inst/include/diversityForestB.h"
#include <RcppEigen.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// divforCpp
Rcpp::List divforCpp(uint treetype, std::string dependent_variable_name, Rcpp::NumericMatrix& input_data, std::vector<std::string> variable_names, uint mtry, uint num_trees, bool verbose, uint seed, uint num_threads, bool write_forest, uint importance_mode_r, uint min_node_size, std::vector<std::vector<double>>& split_select_weights, bool use_split_select_weights, std::vector<std::string>& always_split_variable_names, bool use_always_split_variable_names, std::string status_variable_name, bool prediction_mode, Rcpp::List loaded_forest, Rcpp::RawMatrix snp_data, bool sample_with_replacement, bool probability, std::vector<std::string>& unordered_variable_names, bool use_unordered_variable_names, bool save_memory, uint splitrule_r, std::vector<double>& case_weights, bool use_case_weights, std::vector<double>& class_weights, bool predict_all, bool keep_inbag, std::vector<double>& sample_fraction, double alpha, double minprop, bool holdout, uint prediction_type_r, uint num_random_splits, Eigen::SparseMatrix<double>& sparse_data, bool use_sparse_data, bool order_snps, bool oob_error, uint max_depth, std::vector<std::vector<size_t>>& inbag, bool use_inbag, uint nsplits, uint npairs, double proptry, uint divfortype, std::vector<std::vector<size_t>>& promispairs, uint eim_mode, std::vector<size_t>& metricind);
RcppExport SEXP _diversityForestB_divforCpp(SEXP treetypeSEXP, SEXP dependent_variable_nameSEXP, SEXP input_dataSEXP, SEXP variable_namesSEXP, SEXP mtrySEXP, SEXP num_treesSEXP, SEXP verboseSEXP, SEXP seedSEXP, SEXP num_threadsSEXP, SEXP write_forestSEXP, SEXP importance_mode_rSEXP, SEXP min_node_sizeSEXP, SEXP split_select_weightsSEXP, SEXP use_split_select_weightsSEXP, SEXP always_split_variable_namesSEXP, SEXP use_always_split_variable_namesSEXP, SEXP status_variable_nameSEXP, SEXP prediction_modeSEXP, SEXP loaded_forestSEXP, SEXP snp_dataSEXP, SEXP sample_with_replacementSEXP, SEXP probabilitySEXP, SEXP unordered_variable_namesSEXP, SEXP use_unordered_variable_namesSEXP, SEXP save_memorySEXP, SEXP splitrule_rSEXP, SEXP case_weightsSEXP, SEXP use_case_weightsSEXP, SEXP class_weightsSEXP, SEXP predict_allSEXP, SEXP keep_inbagSEXP, SEXP sample_fractionSEXP, SEXP alphaSEXP, SEXP minpropSEXP, SEXP holdoutSEXP, SEXP prediction_type_rSEXP, SEXP num_random_splitsSEXP, SEXP sparse_dataSEXP, SEXP use_sparse_dataSEXP, SEXP order_snpsSEXP, SEXP oob_errorSEXP, SEXP max_depthSEXP, SEXP inbagSEXP, SEXP use_inbagSEXP, SEXP nsplitsSEXP, SEXP npairsSEXP, SEXP proptrySEXP, SEXP divfortypeSEXP, SEXP promispairsSEXP, SEXP eim_modeSEXP, SEXP metricindSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< uint >::type treetype(treetypeSEXP);
    Rcpp::traits::input_parameter< std::string >::type dependent_variable_name(dependent_variable_nameSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix& >::type input_data(input_dataSEXP);
    Rcpp::traits::input_parameter< std::vector<std::string> >::type variable_names(variable_namesSEXP);
    Rcpp::traits::input_parameter< uint >::type mtry(mtrySEXP);
    Rcpp::traits::input_parameter< uint >::type num_trees(num_treesSEXP);
    Rcpp::traits::input_parameter< bool >::type verbose(verboseSEXP);
    Rcpp::traits::input_parameter< uint >::type seed(seedSEXP);
    Rcpp::traits::input_parameter< uint >::type num_threads(num_threadsSEXP);
    Rcpp::traits::input_parameter< bool >::type write_forest(write_forestSEXP);
    Rcpp::traits::input_parameter< uint >::type importance_mode_r(importance_mode_rSEXP);
    Rcpp::traits::input_parameter< uint >::type min_node_size(min_node_sizeSEXP);
    Rcpp::traits::input_parameter< std::vector<std::vector<double>>& >::type split_select_weights(split_select_weightsSEXP);
    Rcpp::traits::input_parameter< bool >::type use_split_select_weights(use_split_select_weightsSEXP);
    Rcpp::traits::input_parameter< std::vector<std::string>& >::type always_split_variable_names(always_split_variable_namesSEXP);
    Rcpp::traits::input_parameter< bool >::type use_always_split_variable_names(use_always_split_variable_namesSEXP);
    Rcpp::traits::input_parameter< std::string >::type status_variable_name(status_variable_nameSEXP);
    Rcpp::traits::input_parameter< bool >::type prediction_mode(prediction_modeSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type loaded_forest(loaded_forestSEXP);
    Rcpp::traits::input_parameter< Rcpp::RawMatrix >::type snp_data(snp_dataSEXP);
    Rcpp::traits::input_parameter< bool >::type sample_with_replacement(sample_with_replacementSEXP);
    Rcpp::traits::input_parameter< bool >::type probability(probabilitySEXP);
    Rcpp::traits::input_parameter< std::vector<std::string>& >::type unordered_variable_names(unordered_variable_namesSEXP);
    Rcpp::traits::input_parameter< bool >::type use_unordered_variable_names(use_unordered_variable_namesSEXP);
    Rcpp::traits::input_parameter< bool >::type save_memory(save_memorySEXP);
    Rcpp::traits::input_parameter< uint >::type splitrule_r(splitrule_rSEXP);
    Rcpp::traits::input_parameter< std::vector<double>& >::type case_weights(case_weightsSEXP);
    Rcpp::traits::input_parameter< bool >::type use_case_weights(use_case_weightsSEXP);
    Rcpp::traits::input_parameter< std::vector<double>& >::type class_weights(class_weightsSEXP);
    Rcpp::traits::input_parameter< bool >::type predict_all(predict_allSEXP);
    Rcpp::traits::input_parameter< bool >::type keep_inbag(keep_inbagSEXP);
    Rcpp::traits::input_parameter< std::vector<double>& >::type sample_fraction(sample_fractionSEXP);
    Rcpp::traits::input_parameter< double >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< double >::type minprop(minpropSEXP);
    Rcpp::traits::input_parameter< bool >::type holdout(holdoutSEXP);
    Rcpp::traits::input_parameter< uint >::type prediction_type_r(prediction_type_rSEXP);
    Rcpp::traits::input_parameter< uint >::type num_random_splits(num_random_splitsSEXP);
    Rcpp::traits::input_parameter< Eigen::SparseMatrix<double>& >::type sparse_data(sparse_dataSEXP);
    Rcpp::traits::input_parameter< bool >::type use_sparse_data(use_sparse_dataSEXP);
    Rcpp::traits::input_parameter< bool >::type order_snps(order_snpsSEXP);
    Rcpp::traits::input_parameter< bool >::type oob_error(oob_errorSEXP);
    Rcpp::traits::input_parameter< uint >::type max_depth(max_depthSEXP);
    Rcpp::traits::input_parameter< std::vector<std::vector<size_t>>& >::type inbag(inbagSEXP);
    Rcpp::traits::input_parameter< bool >::type use_inbag(use_inbagSEXP);
    Rcpp::traits::input_parameter< uint >::type nsplits(nsplitsSEXP);
    Rcpp::traits::input_parameter< uint >::type npairs(npairsSEXP);
    Rcpp::traits::input_parameter< double >::type proptry(proptrySEXP);
    Rcpp::traits::input_parameter< uint >::type divfortype(divfortypeSEXP);
    Rcpp::traits::input_parameter< std::vector<std::vector<size_t>>& >::type promispairs(promispairsSEXP);
    Rcpp::traits::input_parameter< uint >::type eim_mode(eim_modeSEXP);
    Rcpp::traits::input_parameter< std::vector<size_t>& >::type metricind(metricindSEXP);
    rcpp_result_gen = Rcpp::wrap(divforCpp(treetype, dependent_variable_name, input_data, variable_names, mtry, num_trees, verbose, seed, num_threads, write_forest, importance_mode_r, min_node_size, split_select_weights, use_split_select_weights, always_split_variable_names, use_always_split_variable_names, status_variable_name, prediction_mode, loaded_forest, snp_data, sample_with_replacement, probability, unordered_variable_names, use_unordered_variable_names, save_memory, splitrule_r, case_weights, use_case_weights, class_weights, predict_all, keep_inbag, sample_fraction, alpha, minprop, holdout, prediction_type_r, num_random_splits, sparse_data, use_sparse_data, order_snps, oob_error, max_depth, inbag, use_inbag, nsplits, npairs, proptry, divfortype, promispairs, eim_mode, metricind));
    return rcpp_result_gen;
END_RCPP
}
// numSmaller
Rcpp::IntegerVector numSmaller(Rcpp::NumericVector values, Rcpp::NumericVector reference);
RcppExport SEXP _diversityForestB_numSmaller(SEXP valuesSEXP, SEXP referenceSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type values(valuesSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type reference(referenceSEXP);
    rcpp_result_gen = Rcpp::wrap(numSmaller(values, reference));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_diversityForestB_divforCpp", (DL_FUNC) &_diversityForestB_divforCpp, 51},
    {"_diversityForestB_numSmaller", (DL_FUNC) &_diversityForestB_numSmaller, 2},
    {NULL, NULL, 0}
};

RcppExport void R_init_diversityForestB(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
