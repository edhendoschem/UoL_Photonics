#include <iostream> //Cout
#include <cmath> //Exponential function
#include <fstream> //Create files
#include <algorithm> //Sort algorithm
#include <opencv2/opencv.hpp> //Includes all opencv headers
#include <boost/filesystem.hpp> //Required to iterate over directory

//Libraries required: boost_system;boost_filesystem;opencv_core;opencv_imgproc;opencv_highgui;
//opencv_imgcodecs;

using namespace boost::filesystem;

//Split the name and extension from the path
std::vector<std::string> split_name_ext(std::string const& paths)
{
    std::string const name {paths.substr(0, paths.size()-4)};
    std::string const extension {paths.substr(paths.size()-4, 4)};
    std::vector<std::string> output {name, extension};
    
    return output;
}

//Obtain only the name from the path
std::string get_name(std::string const& paths, char const token = '/')
{
    std::string output;
    for (auto it = std::cend(paths); it != std::cbegin(paths); --it)
    {
        if (*it == token)
        {
            ++it;
            while (it != std::cend(paths))
            {
                output += *it;
                ++it;
            }
            
            break;
        }
    }
    
    output = split_name_ext(output)[0];
    
    return output;
}


int main(int argc, char **argv)
{
    double pixel_length;
    std::cout<<"Please insert pixel length (distance/pixel):\n";
    std::cin>>pixel_length;
    std::cout<<'\n';
    
    path p {current_path()};
    std::vector<std::string> images;
    directory_iterator it{p};
    while (it != directory_iterator{})
    {
        path w = *it;
        std::string st {w.string()};
        //Only add image files
        if (st.substr(st.size()-4, 4) == std::string{".jpg"} ||
            st.substr(st.size()-4, 4) == std::string{".png"} ||
            st.substr(st.size()-4, 4) == std::string{".bmp"})
        {
            images.emplace_back(st);
        }
        ++it;
    }
    
    //Sort by alphabetic order
    std::sort(images.begin(), images.end(),
    [](std::string const& a, std::string const& b)
    {
        std::string a_name {get_name(a)};
        std::string b_name {get_name(b)};
        return a_name < b_name;
    }
    );
    
    
    double pixel_area {pixel_length * pixel_length};
    int loop_count {1};
    
    std::ofstream file_handle {"spot_size.txt"};
    file_handle<<"Pixel length = "<<pixel_length<<'\n';
    file_handle<<"=================================\n";
    
    for (auto i : images)
    {
        std::string location {i};
        //Read the image
        cv::Mat image = cv::imread(location, cv::IMREAD_GRAYSCALE);
        
        //Checking for failure
        if (image.empty())
        {
            std::cout<<"Unable to open image\n";
            return -1;
        }

        //Accept only char type matrices
        CV_Assert(image.depth() == CV_8U);
        const int channels {image.channels()};
        if (channels > 1)
        {
            std::cout<<"Error, multichannel image detected at location:\n"<<i<<'\n';
            return -1;
        }
    
        int rows {image.rows};
        int cols {image.cols};
    
        double max_val {0.0};
        double pixel_count {0.0};
    
        for(auto it = image.begin<uchar>(), end = image.end<uchar>(); it != end; ++it)
        {
            double val {static_cast<double>(*it)};
            max_val = max_val < val ? val : max_val;
        }
    
        //Give margin safety for pixels which may be more intense than background, but not part of 
        //the beam. Use minimum threshold ~5% of max value. max_val must be at least 32
        double threshold {round(max_val * 1.0 / exp(3.0))}; 
    
        for(auto it = image.begin<uchar>(), end = image.end<uchar>(); it != end; ++it)
        {
            double val {static_cast<double>(*it)};
            if (val >= threshold)
            {
                ++pixel_count;
                *it = static_cast<unsigned char> (max_val);
            }
        }
    
        std::vector<std::string> name_ext {split_name_ext(i)};
        std::string name {name_ext[0]+"_"+std::to_string(loop_count)+name_ext[1]};
        
        file_handle<<"Path: "<<i<<'\n';
        file_handle<<"File name: "<<get_name(i)<<'\n';
        file_handle<<"Number of rows = "<<rows<<'\n';
        file_handle<<"Number of columns = "<<cols<<'\n';
        file_handle<<"Number of channels = "<<channels<<'\n';
        file_handle<<"Maximum pixel intensity = "<<max_val<<'\n';
        file_handle<<"Pixel intensity threshold = "<<threshold<<'\n';
        file_handle<<"Pixel count = "<<pixel_count<<'\n';
        file_handle<<"Area = "<<pixel_count * pixel_area<<'\n';
        file_handle<<"___________________________________\n";
        
        bool success {cv::imwrite(name, image)};
    
        if (!success)
        {
            std::cout<<"Error saving image\n";
            return -1;
        }
    
        ++loop_count;
    } //End for loop
    
    file_handle.close();
    

    return 0;
}
