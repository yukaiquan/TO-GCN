#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <algorithm>
#include <cstdlib>  // for exit()

// 显示帮助信息
void show_help(const char* program_name) {
    std::cerr << "Usage: " << program_name << " [options]\n"
              << "Options:\n"
              << "  -i <file>    Input TSV file (required)\n"
              << "  -o <file>    Output file (default: seeds.txt)\n"
              << "  -s <N>       Start column of target condition (0-based, required)\n"
              << "  -e <N>       End column of target condition (0-based, required)\n"
              << "  -m <value>   Min initial expression (default: 10.0)\n"
              << "  -t           Use trend decreasing mode (strict by default)\n"
              << "  -l <slope>   Slope threshold for trend (default: -0.1)\n"
              << "  -h           Show this help\n\n"
              
              << "用法: " << program_name << " [选项]\n"
              << "选项:\n"
              << "  -i <文件>    输入TSV文件路径 (必需)\n"
              << "  -o <文件>    输出文件路径 (默认: seeds.txt)\n"
              << "  -s <N>       目标条件起始列 (从0计数, 必需)\n"
              << "  -e <N>       目标条件结束列 (从0计数, 必需)\n"
              << "  -m <值>      初始时间点最小表达量 (默认: 10.0)\n"
              << "  -t           启用趋势递减模式 (默认严格递减)\n"
              << "  -l <斜率>    趋势判断斜率阈值 (默认: -0.1)\n"
              << "  -h           显示帮助信息\n";
    exit(1);
}
// 判断时间序列是否严格递减
bool is_strictly_decreasing(const std::vector<double>& values) {
    for (size_t i = 1; i < values.size(); ++i) {
        if (values[i] >= values[i-1]) return false;
    }
    return true;
}

// 判断时间序列是否趋势递减
bool is_trend_decreasing(const std::vector<double>& values, double slope_threshold) {
    double sum_x = 0.0, sum_y = 0.0, sum_xy = 0.0, sum_x2 = 0.0;
    const int n = values.size();
    
    for (int i = 0; i < n; ++i) {
        sum_x += i;
        sum_y += values[i];
        sum_xy += i * values[i];
        sum_x2 += i * i;
    }
    
    double slope = (n * sum_xy - sum_x * sum_y) / (n * sum_x2 - sum_x * sum_x);
    return slope < slope_threshold;
}

int main(int argc, char* argv[]) {
    // 默认参数
    std::string input_file, output_file = "seeds.txt";
    int target_start = -1, target_end = -1;
    double min_initial = 10.0, slope_threshold = -0.1;
    bool use_trend = false;

    // 解析命令行参数
    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];
        if (arg == "-i" && i+1 < argc) {
            input_file = argv[++i];
        } else if (arg == "-o" && i+1 < argc) {
            output_file = argv[++i];
        } else if (arg == "-s" && i+1 < argc) {
            target_start = std::stoi(argv[++i]);
        } else if (arg == "-e" && i+1 < argc) {
            target_end = std::stoi(argv[++i]);
        } else if (arg == "-m" && i+1 < argc) {
            min_initial = std::stod(argv[++i]);
        } else if (arg == "-t") {
            use_trend = true;
        } else if (arg == "-l" && i+1 < argc) {
            slope_threshold = std::stod(argv[++i]);
        } else if (arg == "-h") {
            show_help(argv[0]);
        }
    }

    // 验证必需参数
    if (input_file.empty() || target_start < 0 || target_end < target_start) {
        std::cerr << "错误：缺少必需参数或参数无效！\n";
        show_help(argv[0]);
    }

    // 打开文件
    std::ifstream infile(input_file);
    std::ofstream outfile(output_file);
    if (!infile.is_open() || !outfile.is_open()) {
        std::cerr << "错误：无法打开文件！\n";
        return 1;
    }

    // 处理数据
    std::string line;
    while (std::getline(infile, line)) {
        std::istringstream iss(line);
        std::string gene_id;
        std::vector<double> expressions;
        
        if (!(iss >> gene_id)) continue;
        
        double value;
        int col_index = 0;
        while (iss >> value) {
            if (col_index >= target_start && col_index <= target_end) {
                expressions.push_back(value);
            }
            col_index++;
        }
        
        if (expressions.size() != (target_end - target_start + 1)) continue;
        if (expressions[0] < min_initial) continue;
        
        bool is_decreasing = use_trend ? 
            is_trend_decreasing(expressions, slope_threshold) : 
            is_strictly_decreasing(expressions);
        
        if (is_decreasing) {
            outfile << gene_id << std::endl;
        }
    }

    std::cout << "生成种子基因列表: " << output_file << std::endl;
    return 0;
}
