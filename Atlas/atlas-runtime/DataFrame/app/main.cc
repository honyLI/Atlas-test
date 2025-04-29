#include <DataFrame/DataFrame.h>

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <map>
#include <memory>
#include <string>
#include "pointer_shim.h"
#include "pointer.h"
#include "runtime.h"
#include <unistd.h>

using namespace hmdf;

// Download dataset at https://www1.nyc.gov/site/tlc/about/tlc-trip-record-data.page.
// The following code is implemented based on the format of 2016 datasets.
// atlas::AtlasUniquePtr<void> ptrs[1500000][11];

static double haversine(double lat1, double lon1, double lat2, double lon2)
{
    // Distance between latitudes and longitudes
    double dLat = (lat2 - lat1) * M_PI / 180.0;
    double dLon = (lon2 - lon1) * M_PI / 180.0;

    // Convert to radians.
    lat1 = lat1 * M_PI / 180.0;
    lat2 = lat2 * M_PI / 180.0;

    // Apply formulae.
    double a   = pow(sin(dLat / 2), 2) + pow(sin(dLon / 2), 2) * cos(lat1) * cos(lat2);
    double rad = 6371;
    double c   = 2 * asin(sqrt(a));
    return rad * c;
}

// StdDataFrame<uint64_t> load_data()
// {
//     return read_csv<-1, int, SimpleTime, SimpleTime, int, double, double, double, int, char, double,
//                     double, int, double, double, double, double, double, double, double>(
//         "/mnt/all.csv", "vendor_id", "pickup_datetime", "dropoff_datetime",
//         "passenger_count", "trip_distance", "pickup_longitude", "pickup_latitude", "RatecodeID",
//         "store_and_fwd_flag", "dropoff_longitude", "dropoff_latitude", "payment_type",
//         "fare_amount", "extra", "mta_tax", "tip_amount", "tolls_amount", "improvement_surcharge",
//         "total_amount");
// }

StdDataFrame<uint64_t> load_data()
{
    return read_csv<-1,std::string, int, SimpleTime, SimpleTime, int,  double, double, double,
                    double, char, int >(
        "/mnt/all.csv","id", "vendor_id", "pickup_datetime", "dropoff_datetime",
        "passenger_count", "pickup_longitude", "pickup_latitude"
        , "dropoff_longitude", "dropoff_latitude", "store_and_fwd_flag","trip_duration");
}

void print_number_vendor_ids_and_unique(StdDataFrame<uint64_t>& df)
{
    std::cout << "print_number_vendor_ids_and_unique()" << std::endl;
    std::cout << "Number of vendor_ids in the train dataset: " << df.get_column<int>("vendor_id").size() << std::endl;
    // auto vec_ptr = static_cast<std::vector<std::string>*>(ptrs[0][0].deref_get());
    // std::cout << "Number of vendor_ids in the train dataset: " << vec_ptr->size() << std::endl;
    std::cout << "Number of unique vendor_ids in the train dataset:" << df.get_col_unique_values<int>("vendor_id").size() << std::endl;
    std::cout << std::endl;
}

void print_passage_counts_by_vendor_id(StdDataFrame<uint64_t>& df, int vendor_id)
{
    std::cout << "print_passage_counts_by_vendor_id(vendor_id), vendor_id = " << vendor_id << std::endl;

    auto sel_vendor_functor = [&](const uint64_t&, const int& vid) -> bool {
        return vid == vendor_id;
    };
    auto sel_df = df.get_data_by_sel<int, decltype(sel_vendor_functor), int, SimpleTime, double, char>("vendor_id", sel_vendor_functor);
    auto& passage_count_vec = sel_df.get_column<int>("passenger_count");
    std::map<int, int> passage_count_map;
    for (auto passage_count : passage_count_vec) {
        passage_count_map[passage_count]++;
    }
    for (auto& [passage_count, cnt] : passage_count_map) {
        std::cout << "passage_count= " << passage_count << ", cnt = " << cnt << std::endl;
    }
    std::cout << std::endl;
}

void calculate_trip_duration(StdDataFrame<uint64_t>& df)
{
    std::cout << "calculate_trip_duration()" << std::endl;

    auto& pickup_time_vec  = df.get_column<SimpleTime>("pickup_datetime");
    auto& dropoff_time_vec = df.get_column<SimpleTime>("dropoff_datetime");
    assert(pickup_time_vec.size() == dropoff_time_vec.size());

    std::vector<uint64_t> duration_vec;
    duration_vec.reserve(1024*1024*8);
    for (uint64_t i = 0; i < pickup_time_vec.size(); i++) {
        //auto pickup_time_second  = pickup_time_vec[i].to_second();
        SimpleTime * t = static_cast<SimpleTime *>(ptrs[i][2].deref_get());
        auto pickup_time_second = (*t).to_second();
        t = static_cast<SimpleTime *>(ptrs[i][3].deref_get());
        auto dropoff_time_second = (*t).to_second();
        //auto dropoff_time_second = dropoff_time_vec[i].to_second();
        duration_vec.push_back(dropoff_time_second - pickup_time_second);
        ptrs[i][11].reset(&duration_vec.back(),sizeof(uint64_t));
    }

    df.load_column("duration", std::move(duration_vec), nan_policy::dont_pad_with_nans); 
    // MaxVisitor<uint64_t> max_visitor;
    // MinVisitor<uint64_t> min_visitor;
    // MeanVisitor<uint64_t> mean_visitor;
    // df.multi_visit(std::make_pair("duration", &max_visitor),
    //                std::make_pair("duration", &min_visitor),
    //                std::make_pair("duration", &mean_visitor)); //todo
    // std::cout << "Mean duration = " << mean_visitor.get_result() << " seconds" << std::endl;
    // std::cout << "Min duration = " << min_visitor.get_result() << " seconds" << std::endl;
    // std::cout << "Max duration = " << max_visitor.get_result() << " seconds" << std::endl;
    std::cout << std::endl;
}

void calculate_distribution_store_and_fwd_flag(StdDataFrame<uint64_t>& df)
{
    std::cout << "calculate_distribution_store_and_fwd_flag()" << std::endl;

    auto sel_N_saff_functor = [&](const uint64_t&, const char& saff) -> bool {
        return saff == 'N';
    };
    auto N_df =
        df.get_data_by_sel<char, decltype(sel_N_saff_functor), int, SimpleTime, double, char>(
            "store_and_fwd_flag", sel_N_saff_functor);
    std::cout << static_cast<double>(N_df.get_index().size()) / df.get_index().size() << std::endl;

    // auto unique_vendor_id_vec1 = N_df.get_col_unique_values<int>("vendor_id");
    
    // std::cout << '{';
    // for (auto& vector_id : unique_vendor_id_vec1) {
    //     std::cout << vector_id << ", ";
    // }
    // std::cout << '}' << std::endl;

    auto sel_Y_saff_functor = [&](const uint64_t&, const char& saff) -> bool {
        return saff == 'Y';
    };
    auto Y_df =
        df.get_data_by_sel<char, decltype(sel_Y_saff_functor), int, SimpleTime, double, char>(
            "store_and_fwd_flag", sel_Y_saff_functor);
    
    auto unique_vendor_id_vec = Y_df.get_column<int>("vendor_id");
    std::vector<int> result;
    std::unordered_set<int> seen;

    for (int num : unique_vendor_id_vec) {
        if (seen.insert(num).second) {
            result.push_back(num);
        }
    }
    
    std::cout << '{';
    for (auto& vector_id : result) {
        std::cout << vector_id << ", ";
    }
    std::cout << '}' << std::endl;

    std::cout << std::endl;
}

void calculate_haversine_distance_column(StdDataFrame<uint64_t>& df)
{
    std::cout << "calculate_haversine_distance_column()" << std::endl;

    auto& pickup_longitude_vec  = df.get_column<double>("pickup_longitude");
    auto& pickup_latitude_vec   = df.get_column<double>("pickup_latitude");
    auto& dropoff_longitude_vec = df.get_column<double>("dropoff_longitude");
    auto& dropoff_latitude_vec  = df.get_column<double>("dropoff_latitude");
    assert(pickup_longitude_vec.size() == pickup_latitude_vec.size());
    assert(pickup_longitude_vec.size() == dropoff_longitude_vec.size());
    assert(pickup_longitude_vec.size() == dropoff_latitude_vec.size());
    std::vector<double> haversine_distance_vec;
    haversine_distance_vec.reserve(1024*1024*8);
    //std::chrono::time_point<std::chrono::steady_clock> start = std::chrono::steady_clock::now();
    for (uint64_t i = 0; i < pickup_longitude_vec.size(); i++) {

        double pickup_latitude_v = *(static_cast<double *>(ptrs[i][6].deref_get()));
        double pickup_longitude_v = *(static_cast<double *>(ptrs[i][5].deref_get())); 
        double dropoff_latitude_v = *(static_cast<double *>(ptrs[i][8].deref_get()));
        double dropoff_longitude_v = *(static_cast<double *>(ptrs[i][7].deref_get()));

        //double val = haversine(pickup_latitude_vec[i],pickup_longitude_vec[i],dropoff_latitude_vec[i],dropoff_longitude_vec[i]);
        double val = haversine(pickup_latitude_v,pickup_longitude_v,dropoff_latitude_v,dropoff_longitude_v);
        haversine_distance_vec.push_back(val);
        ptrs[i][12].reset(&haversine_distance_vec.back(),sizeof(val));
    }

    //std::chrono::time_point<std::chrono::steady_clock> end = std::chrono::steady_clock::now();
    //std::cout << std::endl;
    //std::cout << "Time taken for haversine calculation: " << std::chrono::duration_cast<std::chrono::microseconds>(end - start).count() << " us" << std::endl;
    //std::cout << std::endl;
    df.load_column("haversine_distance", std::move(haversine_distance_vec), nan_policy::dont_pad_with_nans);

    auto sel_functor = [&](const uint64_t&, const double& dist) -> bool { return dist > 100; };
    auto sel_df = df.get_data_by_sel<double, decltype(sel_functor), int, SimpleTime, double, char>("haversine_distance", sel_functor);
    std::cout << "Number of rows that have haversine_distance > 100 KM = " << sel_df.get_index().size() << std::endl;

    std::cout << std::endl;
}

void analyze_trip_timestamp(StdDataFrame<uint64_t>& df)
{
    std::cout << "analyze_trip_timestamp()" << std::endl;

    MaxVisitor<SimpleTime> max_visitor;
    MinVisitor<SimpleTime> min_visitor;
    df.multi_visit(std::make_pair("pickup_datetime", &max_visitor), std::make_pair("pickup_datetime", &min_visitor));
    std::cout << max_visitor.get_result() << std::endl;
    std::cout << min_visitor.get_result() << std::endl;

    auto& pickup_time_vec = df.get_column<SimpleTime>("pickup_datetime");
    
    std::vector<char> pickup_hour_vec;
    std::vector<char> pickup_day_vec;
    std::vector<char> pickup_month_vec;
    pickup_hour_vec.reserve(1024*1024*16);
    pickup_day_vec.reserve(1024*1024*16);
    pickup_month_vec.reserve(1024*1024*16);
    //std::cout << static_cast<void*>(&(pickup_hour_vec[0])) << std::endl;
    std::map<char, int> pickup_hour_map;
    std::map<char, int> pickup_day_map;
    std::map<char, int> pickup_month_map;

    pickup_hour_vec.resize(pickup_time_vec.size());
    pickup_day_vec.resize(pickup_time_vec.size());
    pickup_month_vec.resize(pickup_time_vec.size());

    auto hour_it  = pickup_hour_vec.begin();
    auto day_it   = pickup_day_vec.begin();
    auto month_it = pickup_month_vec.begin();
    auto time_it  = pickup_time_vec.cbegin();

    for (uint64_t i = 0; i < pickup_time_vec.size();
         ++i, ++hour_it, ++day_it, ++month_it, ++time_it) {
        auto time = *time_it;
        pickup_hour_map[time.hour_]++;
        *hour_it = time.hour_;
        //ptrs[i][13].reset(&(*hour_it),sizeof(time.hour_));
        pickup_day_map[time.day_]++;
        *day_it = time.day_;
        ptrs[i][14].reset(&(*day_it),sizeof(time.day_));
        pickup_month_map[time.month_]++;
        *month_it = time.month_;
        ptrs[i][15].reset(&(*month_it),sizeof(time.month_));
        //std::cout << "time: "<< time.hour_ << " "<< time.day_ << " "<< time.month_ << " " <<std::endl;
    }

    df.load_column("pickup_hour", std::move(pickup_hour_vec), nan_policy::dont_pad_with_nans);
    df.load_column("pickup_day", std::move(pickup_day_vec), nan_policy::dont_pad_with_nans);
    df.load_column("pickup_month", std::move(pickup_month_vec), nan_policy::dont_pad_with_nans);

    //std::cout << "Print top 10 rows." << std::endl;
    //auto top_10_df = df.get_data_by_idx<int, SimpleTime, double, char>(Index2D<StdDataFrame<uint64_t>::IndexType>{0, 9});

    //top_10_df.write<std::ostream, int, SimpleTime, double, char>(std::cout, false, io_format::json);
    //std::cout << std::endl;

    for (auto& [hour, cnt] : pickup_hour_map) {
        std::cout << "pickup_hour = " << static_cast<int>(hour) << ", cnt = " << cnt << std::endl;
    }
    std::cout << std::endl;
    for (auto& [day, cnt] : pickup_day_map) {
        std::cout << "pickup_day = " << static_cast<int>(day) << ", cnt = " << cnt << std::endl;
    }
    std::cout << std::endl;
    for (auto& [month, cnt] : pickup_month_map) {
        std::cout << "pickup_month = " << static_cast<int>(month) << ", cnt = " << cnt << std::endl;
    }
    std::cout << std::endl;
}

// template <typename T_Key>
// void analyze_trip_durations_of_timestamps(StdDataFrame<uint64_t>& df, const char* key_col_name)
// {
//     std::cout << "analyze_trip_durations_of_timestamps() on key = " << key_col_name << std::endl;

//     StdDataFrame<uint64_t> df_key_duration;
//     auto copy_index        = df.get_index();
//     auto copy_key_col      = df.get_column<T_Key>(key_col_name);
//     auto copy_key_duration = df.get_column<uint64_t>("duration");
//     df_key_duration.load_data(std::move(copy_index),
//                               std::make_pair(key_col_name, std::move(copy_key_col)),
//                               std::make_pair("duration", std::move(copy_key_duration)));

//     StdDataFrame<uint64_t> groupby_key =
//         df_key_duration.groupby<GroupbyMedian, T_Key, T_Key, uint64_t>(GroupbyMedian(),key_col_name);

//     auto& key_vec      = groupby_key.get_column<T_Key>(key_col_name);
//     auto& duration_vec = groupby_key.get_column<uint64_t>("duration");
//     for (uint64_t i = 0; i < key_vec.size(); i++) {
//         std::cout << static_cast<int>(key_vec[i]) << " " << duration_vec[i] << std::endl;
//     }

//     std::cout << std::endl;
// }

template <typename T>
T median_aggregator(std::vector<T> vec) {
    if (vec.empty()) return T();
    std::sort(vec.begin(), vec.end());
    size_t mid = vec.size() / 2;
    if (vec.size() % 2 == 0)
        return (vec[mid - 1] + vec[mid]) / T(2);
    else
        return vec[mid];
}

template <typename T_Key>
void analyze_trip_durations_of_timestamps(StdDataFrame<uint64_t>& df, const char* key_col_name)
{
    std::cout << "analyze_trip_durations_of_timestamps() on key = " << key_col_name << std::endl;
    auto it = df.column_tb_.find(key_col_name);
    int col_key = it->second;
    it = df.column_tb_.find("duration");
    int col_dur = it->second;

    auto& key_vec      = df.get_column<T_Key>(key_col_name);     
    auto& duration_vec = df.get_column<uint64_t>("duration");    
    //std::cout <<"key::" << col_key <<" duration::" <<col_dur<< std::endl;
    std::map<T_Key, std::vector<uint64_t>> group_map;

    for (size_t i = 0; i < key_vec.size(); ++i) {   
        T_Key key = *(static_cast<T_Key *>(ptrs[i][col_key].deref_get()));
        //std::cout << 1 << std::endl;
        uint64_t dur = *(static_cast<uint64_t *>(ptrs[i][col_dur].deref_get()));
        group_map[key].push_back(dur);
    }

    for (auto& [key, durations] : group_map) {
        uint64_t median = median_aggregator(durations);
        std::cout << static_cast<int>(key) << " " << median << std::endl;
    }
    std::cout << std::endl;
}

// typedef struct{
//     std::vector<atlas::AtlasUniquePtr<std::vector<std::string>>> id;
//     std::vector<atlas::AtlasUniquePtr<std::vector<int>>> vendor_id;
//     std::vector<atlas::AtlasUniquePtr<std::vector<SimpleTime>>> pickup_datetime;
//     std::vector<atlas::AtlasUniquePtr<std::vector<SimpleTime>>> dropoff_datetime;
//     std::vector<atlas::AtlasUniquePtr<std::vector<int>>> passenger_count;
//     std::vector<atlas::AtlasUniquePtr<std::vector<double>>> pickup_longitude;
//     std::vector<atlas::AtlasUniquePtr<std::vector<double>>> pickup_latitude;
//     std::vector<atlas::AtlasUniquePtr<std::vector<double>>> dropoff_longitude;
//     std::vector<atlas::AtlasUniquePtr<std::vector<double>>> dropoff_latitude;
//     std::vector<atlas::AtlasUniquePtr<std::vector<char>>> store_and_fwd_flag;
//     std::vector<atlas::AtlasUniquePtr<std::vector<int>>> trip_duration;
// } uni_ptr;



int main()
{
    
    runtime_init();

    std::chrono::time_point<std::chrono::steady_clock> times[10];

    auto df  = load_data();   // √ 
    times[0] = std::chrono::steady_clock::now();

    print_number_vendor_ids_and_unique(df);  // √
    times[1] = std::chrono::steady_clock::now();

    print_passage_counts_by_vendor_id(df, 1); // √
    times[2] = std::chrono::steady_clock::now();

    print_passage_counts_by_vendor_id(df, 2); // √
    times[3] = std::chrono::steady_clock::now();

    calculate_trip_duration(df);  //todo
    times[4] = std::chrono::steady_clock::now();

    calculate_distribution_store_and_fwd_flag(df); //√
    times[5] = std::chrono::steady_clock::now();

    calculate_haversine_distance_column(df);  //√
    times[6] = std::chrono::steady_clock::now();

    analyze_trip_timestamp(df);  //√
    times[7] = std::chrono::steady_clock::now();

    analyze_trip_durations_of_timestamps<char>(df, "pickup_day");
    times[8] = std::chrono::steady_clock::now();

    analyze_trip_durations_of_timestamps<char>(df, "pickup_month");
    times[9] = std::chrono::steady_clock::now();
    //sleep(1000);
    for (uint32_t i = 1; i < std::size(times); i++) {
        std::cout << "Step " << i << ": "
                  << std::chrono::duration_cast<std::chrono::microseconds>(times[i] - times[i - 1])
                         .count()
                  << " us" << std::endl;
    }
    std::cout << "Total: "
              << std::chrono::duration_cast<std::chrono::microseconds>(times[9] - times[0]).count()
              << " us" << std::endl;
    runtime_exit();
    return 0;
}
