#pragma once

/* Progressive add smocking patterns,
support adding irragular smocking patch, 
potential app scenario: smocking design on a loose draped cloth
*/
namespace JGSL{

template<int dim = 3>
void progressive_stitching(std::vector<int>& stitchSeq, int framesize, const std::string& filePath = NULL){
    stitchSeq.resize(framesize);
    for(int i =0 ; i < framesize; i++)
        stitchSeq[i] = i;
}

// Progressive adding stitching according to a stithcing index seq. of full smocking pattern 
template<class T, int dim = 3>
void update_stitchingInfo(int cur_frame, const std::vector<int>& stitchSeq, std::vector<VECTOR<int, 3>>& stitchNodes, std::vector<T>& stitchRatio,
    std::vector<VECTOR<int, 3>>& stitchNodes_0, std::vector<T>& stitchRatio_0)
{
    std::cout << "Before stitching size: " << stitchNodes.size() << std::endl;
    std::cout << "Stitching_0 size: " << stitchNodes.size() << std::endl;
    std::cout << "Cur frame: " << cur_frame << std::endl;
    
    stitchNodes.clear();
    stitchRatio.clear();
    // stitchNodes.resize(cur_frame*2);
    // stitchNodes.resize(cur_frame*2);
    std::cout << "After clear size: " << stitchNodes.size() << std::endl;

    stitchNodes.insert(stitchNodes.end(), stitchNodes_0.begin(), stitchNodes_0.begin() + 2*cur_frame);
    stitchRatio.insert(stitchRatio.end(), stitchRatio_0.begin(), stitchRatio_0.begin() + 2*cur_frame);


    // for(int i = 0 ; i < cur_frame; i++){
    //     int stitching_idx = stitchSeq[i];
    //     // stitchNodes.emplace_back(stitchNodes_0[2*stitching_idx]);
    //     // stitchNodes.emplace_back(stitchNodes_0[2*stitching_idx+1]);
    //     // stitchRatio.emplace_back(stitchRatio_0[2*stitching_idx]);
    //     // stitchRatio.emplace_back(stitchRatio_0[2*stitching_idx+1]);
    //     stitchNodes[2*i] = stitchNodes_0[2*i];
    //     stitchNodes[2*i+1] = stitchNodes_0[2*i+1];
    //     stitchRatio[2*i] = stitchRatio_0[2*i];
    //     stitchRatio[2*i+1] = stitchRatio_0[2*i+1];
        
    // }
    std::cout << "Updated stitching size: " << stitchNodes.size() << std::endl;
}

// progressive updating graph


}