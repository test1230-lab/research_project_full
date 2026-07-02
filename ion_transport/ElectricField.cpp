#include "./include/ElectricField.h"
#include <print>

int ElectricField::num_files() const
{
    const auto n = std::ranges::count_if(std::filesystem::directory_iterator(dir),
        [](const auto& entry){ return entry.is_regular_file(); });

    std::print("{} coefficient files in directory\n", n);

    return static_cast<int>(n);
}

//chatgpt code
int ElectricField::extract_number(const std::filesystem::path& p) const
{
    const std::string name = p.stem().string();  // "1-E_100"

    const auto pos = name.find('_');
    if (pos == std::string::npos)
    {
        throw std::runtime_error("Invalid filename: " + name);
    }

    return std::stoi(name.substr(pos + 1));
}

void ElectricField::read_coeffs()
{
    std::vector<std::pair<int, std::filesystem::path>> files;

    for (const auto& entry : std::filesystem::directory_iterator(dir))
    {
        if (entry.is_regular_file())
        {
            const int num = extract_number(entry.path());
            e_field_vals.push_back(static_cast<double>(num));
            files.push_back({num, entry.path()});
        }
    }

    std::ranges::sort(files);
    std::ranges::sort(e_field_vals);

    int n = 0;
    for (const auto& file : files)
    {
        std::filesystem::directory_entry entry(file.second);
        if (entry.is_regular_file()) 
        {
            std::ifstream in(entry.path());
            std::string line;

            //skip line 1,2
            std::getline(in, line);
            std::getline(in, line);

            //get ion thermal speeds
            std::getline(in, line);
            int idx = 0;
            for (auto&& elem : std::views::split(line, delim))
            {
                const double val = std::stod(std::string(elem.begin(), elem.end()));
                ion_thermal_speeds[idx++, n] = val;
            }

            //get coeffs
            for (int i = 0; i < n_cols; i++)
            {
                std::getline(in, line);
                idx = 0;
                for (auto&& elem : std::views::split(line, delim))
                {
                    const double val = std::stod(std::string(elem.begin(), elem.end()));
                    coeffs[idx++, n, i] = val;
                }
            }

            n++;
        }
    }
}

void ElectricField::read_coeffs_new_fmt()
{
    std::vector<std::pair<int, std::filesystem::path>> files;

    for (const auto& entry : std::filesystem::directory_iterator(dir))
    {
        if (entry.is_regular_file())
        {
            const int num = extract_number(entry.path());
            e_field_vals.push_back(static_cast<double>(num));
            files.push_back({num, entry.path()});
        }
    }

    std::ranges::sort(files);
    std::ranges::sort(e_field_vals);

    int n = 0;
    for (const auto& file : files)
    {
        std::filesystem::directory_entry entry(file.second);
        if (entry.is_regular_file()) 
        {
            std::ifstream in(entry.path());
            std::string line;

            //skip line 1-8 inclusive
            for (int i = 0; i < 8; i++)
            {
                std::getline(in, line);
            }

            //get ion thermal speeds
            std::getline(in, line);
            int idx = 0;
            for (auto&& elem : std::views::split(line, delim))
            {
                const double val = std::stod(std::string(elem.begin(), elem.end()));
                ion_thermal_speeds[idx++, n] = val;
            }

            //get coeffs
            for (int i = 0; i < 2*n_cols - 1; i++)
            {
                std::getline(in, line);
                if (i % 2 == 0)
                {
                    const int coeff_idx = i / 2;
                    idx = 0;
                    for (auto&& elem : std::views::split(line, delim))
                    {
                        const double val = std::stod(std::string(elem.begin(), elem.end()));
                        coeffs[idx++, n, coeff_idx] = val;
                        
                    }
                }
            }
            

            n++;
        }
    }
}

double ElectricField::compute_dist_discrete(int electric_field, int aspect_angle, double v) const
{
    const int aspect_angle_idx = aspect_angle / 10;

    auto it = std::ranges::lower_bound(e_field_vals, static_cast<double>(electric_field));

    if (it == e_field_vals.end() || *it != static_cast<double>(electric_field))
    {
        throw std::runtime_error("Electric field not found in coefficient table");
    }

    const int e_idx = static_cast<int>(it - e_field_vals.begin());

    //y is vx/b, where vx is the line-of-sight speed and b is the ion thermal speed).
    const double y = v/ion_thermal_speeds[aspect_angle_idx, e_idx];
    const double x = y/4.0;

    if (std::abs(x) > 1.0)
    {
        return 0.0;
    }

    double sum = 0.0;
    for (int i = 0; i < n_cols; i++)
    {
        const double ci = coeffs[aspect_angle_idx, e_idx, i];
        sum += ci*std::legendre(i*2, x);
    }

    return sum/ion_thermal_speeds[aspect_angle_idx, e_idx];

    
}

double ElectricField::compute_dist(double electric_field, int aspect_angle, double v) const
{
    const int angle_idx = aspect_angle / 10;
    electric_field = std::clamp(electric_field, 0.0, 200.0);
    const double e_norm = 2.0*((electric_field - e_field_vals[0])/(e_field_vals[n_files - 1] - e_field_vals[0])) - 1.0;

    const double ion_thermal_speed_interp = (*ion_speed_interp[angle_idx])(e_norm);
    const double y = v/ion_thermal_speed_interp;
    const double x = y/4.0;

    if (std::abs(x) > 1.0)
    {
        return 0.0;
    }

    std::array<double, n_cols> coeffs1;
    for (int i = 0; i < n_cols; i++)
    {
        coeffs1[i] = (*e_interp[angle_idx][i])(e_norm);
    }

    return eval_legendre_series(coeffs1, x)/ion_thermal_speed_interp;
}

double ElectricField::eval_legendre_series(const std::array<double, n_cols>& c, double x) const 
{
    double sum = c[0];
    double pm1 = 1.0;   // P_{i-1}
    double p   = x;     // P_i  (i = 1)

    for (int i = 1; i < 2*(n_cols - 1); i++) 
    {
        const double pnext = leg_a[i]*x*p - leg_b[i]*pm1;  //no division
        pm1 = p; p = pnext;
        if (i % 2 == 1)
        {
            sum = std::fma(p, c[(i + 1)/2], sum);
        }
    }

    return sum;
}

void ElectricField::compute_interp_coeffs()
{
    std::vector<double> e_norm(n_files);

    for (int i = 0; i < n_files; i++)
    {
        e_norm[i] = 2.0*((e_field_vals[i] - e_field_vals[0])/(e_field_vals[n_files - 1] - e_field_vals[0])) - 1.0;
    }

    for (int i = 0; i < n_angles; i++)
    {
        for (int j = 0; j < n_cols; j++)
        {
            std::vector<double> x = e_norm;
            std::vector<double> y(n_files);

            for (int k = 0; k < n_files; k++)
            {
                y[k] = coeffs[i, k, j];
            }

            e_interp[i][j].emplace(makima(std::move(x), std::move(y)));
        }   
    }

    for (int i = 0; i < n_angles; i++)
    {
        std::vector<double> x = e_norm;
        std::vector<double> y(n_files);

        for (int j = 0; j < n_files; j++)
        {
            y[j] = ion_thermal_speeds[i, j];
        }

        ion_speed_interp[i].emplace(pchip(std::move(x), std::move(y)));
    }

}

double ElectricField::compute_integral(double electric_field, int aspect_angle) const
{
    const int angle_idx = aspect_angle / 10;
    electric_field = std::clamp(electric_field, 0.0, 200.0);
     const double e_norm = 2.0*((electric_field - e_field_vals[0])/(e_field_vals[n_files - 1] - e_field_vals[0])) - 1.0;
    return 8.0*(*e_interp[angle_idx][0])(e_norm);
}