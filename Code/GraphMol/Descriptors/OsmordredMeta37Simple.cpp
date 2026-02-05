// Simple Meta37 Integration - Direct Osmordred + MW features
// Uses feature selection for BP, MP, Flashpoint (Fixed models)

#include "Osmordred.h"
#include "OsmordredMeta37Dd.h"
#include "OsmordredMeta37Dh.h"
#include "OsmordredMeta37Dp.h"
#include "OsmordredMeta37BpFixed.h"       // Uses feature selection
#include "OsmordredMeta37Logvp.h"
#include "OsmordredMeta37Logpow.h"
#include "OsmordredMeta37Logws.h"
#include "OsmordredMeta37Deltahf.h"
#include "OsmordredMeta37Deltahc.h"
#include "OsmordredMeta37MpFixed.h"       // Uses feature selection
#include "OsmordredMeta37FlashpointFixed.h"  // Uses feature selection
#include "OsmordredMeta37Loghenrycc.h"
#include "OsmordredMeta37Dipolemoment.h"

#include <GraphMol/ROMol.h>
#include <GraphMol/Descriptors/MolDescriptors.h>
#include <vector>
#include <cmath>
#include <fstream>
#include <sstream>
#include <chrono>
#include <thread>
#include <mutex>

// #region agent log - Debug instrumentation for cascade bug investigation
namespace {
    // Thread-safe logging helper
    inline void logDebug(const std::string& location, const std::string& message, 
                        const std::string& data, const std::string& hypothesisId) {
        static std::mutex logMutex;
        std::lock_guard<std::mutex> lock(logMutex);
        
        std::ofstream logFile("/Users/guillaume-osmo/Github/osmomain/.cursor/debug.log", 
                             std::ios::app);
        if (logFile.is_open()) {
            auto now = std::chrono::system_clock::now();
            auto timestamp = std::chrono::duration_cast<std::chrono::milliseconds>(
                now.time_since_epoch()).count();
            
            std::stringstream json;
            json << "{\"id\":\"log_" << timestamp << "_" << std::this_thread::get_id()
                 << "\",\"timestamp\":" << timestamp
                 << ",\"location\":\"" << location << "\""
                 << ",\"message\":\"" << message << "\"";
            if (!data.empty()) {
                json << ",\"data\":" << data;
            }
            json << ",\"sessionId\":\"debug-session\""
                 << ",\"runId\":\"run1\"";
            if (!hypothesisId.empty()) {
                json << ",\"hypothesisId\":\"" << hypothesisId << "\"";
            }
            json << "}\n";
            
            logFile << json.str();
            logFile.close();
        }
    }
}
// #endregion agent log

namespace RDKit {
namespace Descriptors {
namespace Osmordred {

// Build features: Osmordred (3585) + MW (1) = 3586 total
std::vector<double> buildSimpleFeatures(const ROMol& mol) {
    // #region agent log
    logDebug("OsmordredMeta37Simple.cpp:buildSimpleFeatures", "Function entry", 
             "{\"thread_id\":" + std::to_string(std::hash<std::thread::id>{}(std::this_thread::get_id())) + "}",
             "C");
    // #endregion agent log
    
    std::vector<double> osmo = calcOsmordred(mol);
    
    // #region agent log
    std::stringstream osmoStr;
    osmoStr << "{\"osmo_size\":" << osmo.size();
    if (osmo.size() > 0) {
        osmoStr << ",\"osmo_sample\":[" << osmo[0];
        for (size_t i = 1; i < std::min(osmo.size(), size_t(5)); ++i) {
            osmoStr << "," << osmo[i];
        }
        osmoStr << "]";
    }
    osmoStr << "}";
    logDebug("OsmordredMeta37Simple.cpp:buildSimpleFeatures", "After calcOsmordred", 
             osmoStr.str(), "C");
    // #endregion agent log
    
    double MW = calcAMW(mol);
    
    std::vector<double> features;
    features.reserve(3586);
    
    // Copy Osmordred
    for (double v : osmo) {
        if (std::isinf(v)) v = 0.0;
        if (std::isnan(v)) v = 0.0;
        features.push_back(v);
    }
    
    // Add MW
    features.push_back(MW);
    
    // #region agent log
    logDebug("OsmordredMeta37Simple.cpp:buildSimpleFeatures", "Function exit", 
             "{\"features_size\":" + std::to_string(features.size()) + 
             ",\"MW\":" + std::to_string(MW) + "}", "C");
    // #endregion agent log
    
    return features;
}

// Calculate all Meta37 properties using simple features
std::vector<double> calcMeta37Simple(const ROMol& mol) {
    // #region agent log
    logDebug("OsmordredMeta37Simple.cpp:calcMeta37Simple", "Function entry", 
             "{\"thread_id\":" + std::to_string(std::hash<std::thread::id>{}(std::this_thread::get_id())) + "}",
             "C");
    // #endregion agent log
    
    std::vector<double> features = buildSimpleFeatures(mol);
    const double* f = features.data();
    std::vector<double> results(13);
    
    // #region agent log
    logDebug("OsmordredMeta37Simple.cpp:calcMeta37Simple", "Before predictions", 
             "{\"features_size\":" + std::to_string(features.size()) + "}", "C");
    // #endregion agent log
    
    results[0] = ::Osmordred::Meta37Dd::predict(f);           // dD
    results[1] = ::Osmordred::Meta37Dh::predict(f);           // dH
    results[2] = ::Osmordred::Meta37Dp::predict(f);           // dP
    results[3] = ::Osmordred::Meta37BpFixed::predictFromFull(features);    // BP (feature selection)
    results[4] = ::Osmordred::Meta37Logvp::predict(f);        // logVP
    results[5] = ::Osmordred::Meta37Logpow::predict(f);       // logPow
    
    // #region agent log
    logDebug("OsmordredMeta37Simple.cpp:calcMeta37Simple", "After first 6 predictions", 
             "{\"dD\":" + std::to_string(results[0]) + ",\"dH\":" + std::to_string(results[1]) + 
             ",\"dP\":" + std::to_string(results[2]) + ",\"BP\":" + std::to_string(results[3]) + 
             ",\"logVP\":" + std::to_string(results[4]) + ",\"logPow\":" + std::to_string(results[5]) + "}", 
             "C");
    // #endregion agent log
    
    results[6] = ::Osmordred::Meta37Logws::predict(f);        // logWS
    results[7] = ::Osmordred::Meta37Deltahf::predict(f);      // deltaHf
    results[8] = ::Osmordred::Meta37Deltahc::predict(f);      // deltaHc
    results[9] = ::Osmordred::Meta37MpFixed::predictFromFull(features);    // MP (feature selection)
    results[10] = ::Osmordred::Meta37FlashpointFixed::predictFromFull(features);  // FP (feature selection)
    results[11] = ::Osmordred::Meta37Loghenrycc::predict(f);  // logHenrycc
    results[12] = ::Osmordred::Meta37Dipolemoment::predict(f);// dipolemoment
    
    // #region agent log
    std::stringstream resultsStr;
    resultsStr << "[";
    for (size_t i = 0; i < results.size(); ++i) {
        if (i > 0) resultsStr << ",";
        resultsStr << results[i];
    }
    resultsStr << "]";
    logDebug("OsmordredMeta37Simple.cpp:calcMeta37Simple", "Function exit", 
             "{\"results\":" + resultsStr.str() + "}", "C");
    // #endregion agent log
    
    return results;
}

}  // namespace Osmordred
}  // namespace Descriptors
}  // namespace RDKit
