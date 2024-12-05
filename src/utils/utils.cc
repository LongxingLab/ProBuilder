#include "utils/utils.hh"
#include "basic/types.hh"
#include "scene/Pose.hh"
#include "sampling/FragMover.hh"
#include "utils/random_util.hh"

#include <cmath>
#include <ctime>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <random>
#include <algorithm>

namespace utils {

using namespace basic;

bool check_break(scene::Pose const & pose)
{
    Size len = pose.size();

    bool flag = false;
    for(Size ires=1; ires<len; ++ires)
    {
        Vec atom_C = pose.xyz(ires, ATOM_C);
        Vec atom_N = pose.xyz(ires+1, ATOM_N);

        Real dist = (atom_C - atom_N).squaredNorm();
        
        if(dist > 2.5) {
            flag = true;
            break;
        }
    }
    return flag;
}
std::pair<bool,std::string> f1 (Size res_len,std::string ss_now,char this_ss,Size min_helix_len,Size max_helix_len,Size min_loop_len,Size max_loop_len,Size now_segment=1,Size res_segment=-1){
        std::uniform_int_distribution<Size> loop_len_random( min_loop_len,max_loop_len);
        // let the first helix and the last helix longer than 16
        // std::uniform_int_distribution<Size> helix_len_random(now_segment!=1?min_helix_len:ceil((min_helix_len+max_helix_len)/2.0),max_helix_len);
        std::uniform_int_distribution<Size> helix_len_random(now_segment!=1?min_helix_len:16,max_helix_len);
        char next_ss = this_ss=='H'?'L':'H';
        std::string ss_append;
        assert(res_segment>=-1);
        ss_append = (this_ss=='H'?std::string(helix_len_random(rng()),'H'):std::string(loop_len_random(rng()),'L'));
        // if(res_segment==0)return std::pair<bool,std::string>(true,ss_now);
        // if(res_segment!=-1 && res_segment>0) return f1(res_len-ss_append.size(),ss_now+ss_append,next_ss,min_helix_len,max_helix_len,min_loop_len, max_loop_len,res_segment-1);
        
        if(res_len<=0)return std::pair<bool,std::string>(false,"H");
        if(this_ss=='H'&& res_len<min_helix_len)return std::pair<bool,std::string>(false,"H");
        if(this_ss=='H'&& res_len<max_helix_len){
            std::uniform_int_distribution<Size> tolerant_len(ceil((min_helix_len+max_helix_len)/2.0),max_helix_len);
            if(res_len>=tolerant_len(rng()))return std::pair<bool,std::string>(true,ss_now+std::string(res_len,'H'));
            return std::pair<bool,std::string>(false,"H");
        }
        if(this_ss=='L'&& res_len<min_loop_len) return std::pair<bool,std::string>(false,"H");
        return f1(res_len-ss_append.size(),ss_now+ss_append,next_ss,min_helix_len,max_helix_len,min_loop_len, max_loop_len,now_segment+1,res_segment>0?(res_segment-1):res_segment);
}

std::string random_helix_dssp(Size min_helix_len,Size max_helix_len,Size min_loop_len,Size max_loop_len,Size max_len,Size segment_num){
    std::string ss;
    std::pair<bool,std::string> result;
    int count=0;
    while (true){
        count++;
        if(segment_num!=-1)result = f1(-1,"",'H',min_helix_len,max_helix_len,min_loop_len,max_loop_len,segment_num);
        else result = f1(max_len,"",'H',min_helix_len,max_helix_len,min_loop_len,max_loop_len);
        if(result.first)return result.second;
        if(count>10000){
            std::cout<<"might some bug in dssp generation,check option is reasonable"<<std::endl;
        }
    }  
}

bool random_bundle_dssp(Size length, Size num_helix, Size min_helix_len, Size max_helix_len, Size min_loop_len, Size max_loop_len, std::string & dssp_prefix, std::string & dssp)
{

    bool pass = true;
    std::vector<Size> loops(num_helix-1, 0);
    std::vector<Size> helices(num_helix, 0);

    Size max_tries = 100000;

    Size min_terminal_helix_len = int( (length - (num_helix-1)*4)/num_helix );
    if (min_terminal_helix_len>14) min_terminal_helix_len = 14;

    while(max_tries--) {
        Size res_len = length;
        
        for(Size idx=0; idx<loops.size(); ++idx) {
            Size this_loop_len = random_int(min_loop_len, max_loop_len);
            loops[idx] = this_loop_len;
            res_len -= this_loop_len;
        }

        for(Size idx=0; idx<helices.size(); ++idx) helices[idx] = 0;
        while(res_len--) {
            Size idx = random_int(0, num_helix-1);
            helices[idx] += 1;
        }

        pass = true;
        if (helices[0] < min_terminal_helix_len || helices[num_helix-1] < min_terminal_helix_len) pass = false;
        for(Size idx=0; idx<helices.size(); ++idx) { if(helices[idx]<min_helix_len || helices[idx]>max_helix_len) pass = false; }
        
        if(pass) break;
    }

    dssp_prefix = "";
    dssp = "";
    std::stringstream prefix;

    for(Size idx=0; idx<num_helix-1; idx++) {
        dssp += std::string(helices[idx], 'H') + std::string(loops[idx], 'L');
        prefix << "H" << helices[idx] << "L" << loops[idx];
    }
    dssp += std::string(helices[num_helix-1], 'H');
    prefix << "H" << helices[num_helix-1];

    dssp_prefix = prefix.str();

    return pass;

}

bool random_bundle_dssp_AzoF(Size length, Size num_helix, Size min_helix_len, Size max_helix_len, Size min_loop_len, Size max_loop_len, std::string & dssp_prefix, std::string & dssp, Size & insert_pos)
{

    bool pass = true;
    std::vector<Size> loops(num_helix-1, 0);
    std::vector<Size> helices(num_helix, 0);

    Size max_tries = 100000;

    Size min_terminal_helix_len = int( (length - (num_helix-1)*4)/num_helix );
    if (min_terminal_helix_len>14) min_terminal_helix_len = 14;

    while(max_tries--) {
        Size res_len = length;
        
        for(Size idx=0; idx<loops.size(); ++idx) {
            Size this_loop_len = random_int(min_loop_len, max_loop_len);
            loops[idx] = this_loop_len;
            res_len -= this_loop_len;
        }

        for(Size idx=0; idx<helices.size(); ++idx) helices[idx] = 0;
        while(res_len--) {
            Size idx = random_int(0, num_helix-1);
            helices[idx] += 1;
        }

        pass = true;
        if (helices[0] < min_terminal_helix_len || helices[num_helix-1] < min_terminal_helix_len) pass = false;
        for(Size idx=0; idx<helices.size(); ++idx) { if(helices[idx]<min_helix_len || helices[idx]>max_helix_len) pass = false; }

        if( random_real(0.0, 1.0) > 0.5 && helices[0] > 20 ) {
            // the first helix
            insert_pos = random_int( helices[0]-5, helices[0] );
        } else if (helices[num_helix-1]>20) {
            // the last helix
            insert_pos = random_int(length-helices[num_helix-1]+1+1,length-helices[num_helix-1]+1+7);
        } else {
            pass = false;
        }
        
        if(pass) break;
    }

    dssp_prefix = "";
    dssp = "";
    std::stringstream prefix;

    for(Size idx=0; idx<num_helix-1; idx++) {
        dssp += std::string(helices[idx], 'H') + std::string(loops[idx], 'L');
        prefix << "H" << helices[idx] << "L" << loops[idx];
    }
    dssp += std::string(helices[num_helix-1], 'H');
    prefix << "H" << helices[num_helix-1];

    dssp_prefix = prefix.str();

    return pass;

}

bool random_bundle_dssp_HQA(Size length, Size num_helix, Size min_helix_len, Size max_helix_len, Size min_loop_len, Size max_loop_len, std::string & dssp_prefix, std::string & dssp, Size & insert_pos)
{

    bool pass = true;
    std::vector<Size> loops(num_helix-1, 0);
    std::vector<Size> helices(num_helix, 0);

    Size max_tries = 100000;

    // hard code everything here
    min_helix_len = 12;

    Size min_terminal_helix_len = int( (length - (num_helix-1)*4)/num_helix );
    if (min_terminal_helix_len>14) min_terminal_helix_len = 14;

    while(max_tries--) {
        Size res_len = length;
        
        for(Size idx=0; idx<loops.size(); ++idx) {
            Size this_loop_len = random_int(min_loop_len, max_loop_len);
            loops[idx] = this_loop_len;
            res_len -= this_loop_len;
        }

        for(Size idx=0; idx<helices.size(); ++idx) helices[idx] = 0;
        while(res_len--) {
            Size idx = random_int(0, num_helix-1);
            helices[idx] += 1;
        }

        pass = true;
        if (helices[0] < min_terminal_helix_len || helices[num_helix-1] < min_terminal_helix_len) pass = false;
        for(Size idx=0; idx<helices.size(); ++idx) { if(helices[idx]<min_helix_len || helices[idx]>max_helix_len) pass = false; }
        
        if(pass) break;
    }


    // pick the insertion pos
    Size pick_helix = random_int(0, num_helix-1);
    Size pick_position = random_int(1, helices[pick_helix]-10);
    insert_pos = 0;
    for(Size idx=0; idx<pick_helix; ++idx) {
        insert_pos += helices[idx] + loops[idx];
    }
    insert_pos += pick_position + 5;

    dssp_prefix = "";
    dssp = "";
    std::stringstream prefix;

    for(Size idx=0; idx<num_helix-1; idx++) {
        dssp += std::string(helices[idx], 'H') + std::string(loops[idx], 'L');
        prefix << "H" << helices[idx] << "L" << loops[idx];
    }
    dssp += std::string(helices[num_helix-1], 'H');
    prefix << "H" << helices[num_helix-1];

    dssp_prefix = prefix.str();

    return pass;

}

bool random_motif_bundle_dssp(Size length,std::string motif_ss,bool side_require, Size num_helix, Size min_helix_len, Size max_helix_len, Size min_loop_len, Size max_loop_len, std::string & dssp_prefix, std::string & dssp,Size & insert_pos)
{

    bool pass = true;
    std::vector<Size> loops(num_helix-1, 0);
    std::vector<Size> helices(num_helix, 0);
    if(side_require){
        motif_ss.insert(motif_ss.begin(),motif_ss[0]);
        motif_ss.push_back(motif_ss.back());
    }
    Size max_tries = 100000;

    // Size min_terminal_helix_len = int( (length - (num_helix-1)*4)/num_helix );
    // if (min_terminal_helix_len>14) min_terminal_helix_len = 14;

    while(max_tries--) {
        Size res_len = length;
        
        for(Size idx=0; idx<loops.size(); ++idx) {
            Size this_loop_len = random_int(min_loop_len, max_loop_len);
            loops[idx] = this_loop_len;
            res_len -= this_loop_len;
        }

        for(Size idx=0; idx<helices.size(); ++idx) helices[idx] = 0;
        while(res_len--) {
            Size idx = random_int(0, num_helix-1);
            helices[idx] += 1;
        }

        pass = true;
        if (helices[0] < 14 || helices[num_helix-1] < 14) pass = false;
        for(Size idx=0; idx<helices.size(); ++idx) { if(helices[idx]<min_helix_len || helices[idx]>max_helix_len) pass = false; }
        dssp = "";
        for(Size idx=0; idx<num_helix-1; idx++) {
            dssp += std::string(helices[idx], 'H') + std::string(loops[idx], 'L');
        }
        dssp += std::string(helices[num_helix-1], 'H');
        if(dssp.find(motif_ss)==std::string::npos)pass=false;
        if(pass) break;
    }


    dssp_prefix = "";
    dssp = "";
    std::stringstream prefix;

    for(Size idx=0; idx<num_helix-1; idx++) {
        dssp += std::string(helices[idx], 'H') + std::string(loops[idx], 'L');
        prefix << "H" << helices[idx] << "L" << loops[idx];
    }
    dssp += std::string(helices[num_helix-1], 'H');
    prefix << "H" << helices[num_helix-1];

    dssp_prefix = prefix.str();

    // find possible insert position
    std::vector<Size> candidate_pos;
    for(Size ipos=0; ipos < length-motif_ss.length(); ++ipos) {
        if(dssp.substr(ipos, motif_ss.length()) == motif_ss) {

            if(side_require) {
                candidate_pos.push_back(ipos+2);
            } else {
                candidate_pos.push_back(ipos+1);
            }

        }
    }
    //for(std::size_t i =0;i!=std::string::npos;){
    //    i = (i==0?i:i+1);
    //    i=dssp.find(motif_ss,i);
    //    if(i!=std::string::npos)candidate_pos.push_back(side_require?i+2:i+1);
    //}
    insert_pos = candidate_pos[rand()%candidate_pos.size()];
    return pass;

}

bool random_motif_pair_bundle_dssp(Size length,std::string motif1_ss, std::string motif2_ss, Size num_helix, Size min_helix_len, Size max_helix_len, Size min_loop_len, Size max_loop_len, std::string & dssp_prefix, std::string & dssp, std::string & fake_dssp, Size & insert_pos1, Size & insert_pos2)
{

    const Size ter_min_helix_len = 14;
    const Size motif_min_sep = 1;

    bool pass = false;
    Size motif1_insert_pos(-1);
    Size motif2_insert_pos(-1);
    // std::string ss;

    const Size motif1_len = motif1_ss.length();
    const Size motif2_len = motif2_ss.length();

    Size max_tries = 100000;
    while(max_tries--) {
        // initialize the string
        std::string work_ss(length, 'X');

        // place motif #1
        Size motif1_first_L = motif1_ss.find_first_of('L');
        Size motif1_last_L = motif1_ss.find_last_of('L');

        if(motif1_first_L == std::string::npos) {
            // all helix 
            motif1_insert_pos = random_int(0, length-motif1_len);
        } else {
            // has loop motif
            Size start_pos = std::max(Size(0), ter_min_helix_len-motif1_first_L);
            Size stop_pos = length-motif1_len - std::max(Size(0), ter_min_helix_len-(motif1_len-(motif1_last_L+1)));
            motif1_insert_pos = random_int(start_pos, stop_pos);
        }
        for(Size idx=0; idx<motif1_len; ++idx) {
            work_ss[motif1_insert_pos+idx] = motif1_ss[idx];
        }

        // now place the motif #2
        std::vector<Size> motif2_allowed_pos;
        Size motif2_first_L = motif2_ss.find_first_of('L');
        Size motif2_last_L = motif2_ss.find_last_of('L');
        if(motif2_first_L == std::string::npos) {
            if(motif1_insert_pos-motif2_len>=motif_min_sep) {
                for(Size idx=0; idx<=motif1_insert_pos-motif2_len-motif_min_sep; ++idx) {
                    motif2_allowed_pos.push_back(idx);
                }
            }
            if(motif1_insert_pos+motif1_len+motif_min_sep+motif2_len<=length) {
                for(Size idx=motif1_insert_pos+motif1_len+motif_min_sep; idx<=length-motif2_len; ++idx) {
                    motif2_allowed_pos.push_back(idx);
                }
            }
        } else {
            if(motif1_insert_pos-motif2_len>=motif_min_sep) {
                for(Size idx=std::max(Size(0), ter_min_helix_len-motif2_first_L); idx<=motif1_insert_pos-motif2_len-motif_min_sep; ++idx) {
                    motif2_allowed_pos.push_back(idx);
                }
            }
            if(motif1_insert_pos+motif1_len+motif_min_sep+motif2_len<=length) {
                for(Size idx=motif1_insert_pos+motif1_len+motif_min_sep; idx<=length-motif2_len- std::max(Size(0), ter_min_helix_len-(motif2_len-(motif2_last_L+1))); ++idx) {
                    motif2_allowed_pos.push_back(idx);
                }
            }
        }

        // finally chose the insertion pos of motif2 and replace the ss str
        if(motif2_allowed_pos.size() == 0) {
            continue;
        } else {
            motif2_insert_pos = motif2_allowed_pos[rand()%motif2_allowed_pos.size()];
        }
        for(Size idx=0; idx<motif2_len; ++idx) {
            work_ss[motif2_insert_pos+idx] = motif2_ss[idx];
        }

        Size cur_num_loops(0);
        for(Size idx=1; idx<length; ++idx) {
            if(work_ss[idx]!='L' && work_ss[idx-1]=='L') {
                cur_num_loops +=1;
            }
        }

        // global check flag
        bool err = false;

        // now add loops
        for(Size iloop=0; iloop<num_helix-1-cur_num_loops; ++iloop) {

            Size this_loop_len = random_int(min_loop_len, max_loop_len);
            std::vector<Size> loop_allowed_pos;
            for(Size idx=ter_min_helix_len; idx<=length-this_loop_len-ter_min_helix_len; ++idx) {
                bool flag = true;
                // overlap with others?
                for(Size jdx=0; jdx<this_loop_len; ++jdx) {
                    if(work_ss[idx+jdx] != 'X') {
                        flag = false;
                        break;
                    }
                }
                // check left
                for(Size jdx=idx-1; jdx>=0 && jdx>= idx-min_helix_len; --jdx) {
                    if(work_ss[jdx] == 'L') {
                        flag = false;
                        break;
                    }
                }
                // check right
                for(Size jdx=idx+1; jdx<length && jdx< idx+this_loop_len+min_helix_len; ++jdx) {
                    if(work_ss[jdx] == 'L') {
                        flag = false;
                        break;
                    }
                }
                if(flag) {
                    loop_allowed_pos.push_back(idx);
                }
            }
            if(loop_allowed_pos.size() == 0) {
                err = true;
                break;
            } else {
                Size this_loop_insert_pos = loop_allowed_pos[rand()%loop_allowed_pos.size()];
                for(Size idx=0; idx<this_loop_len; ++idx) {
                    work_ss[this_loop_insert_pos+idx] = 'L';
                }
            }
        }
        if(err) {
            continue;
        }


        // check the helix length
        // max_helix_len check
        //
        Size temp_helix_len(0);
        for(Size idx=0; idx<length; ++idx) {
            if(work_ss[idx] != 'L') {
                temp_helix_len += 1;
            } else {
                if(temp_helix_len > max_helix_len) {
                    err = true;
                }
                temp_helix_len = 0;
            }
        }
        if(err) {
            continue;
        }

        dssp = work_ss;
        pass = true;
        break;
    }


    // get real dssp
    // fake dssp
    // dssp prefix
    insert_pos1 = motif1_insert_pos+1; // 1 index
    insert_pos2 = motif2_insert_pos+1; // 1 index

    for(Size idx=0; idx<length; ++idx) {
        if(dssp[idx] == 'X') {
            dssp[idx] = 'H';
        }
    }

    fake_dssp = dssp;
    for(Size idx=0; idx<motif1_len; ++idx) {
        fake_dssp[idx+motif1_insert_pos] = 'X';
    }
    for(Size idx=0; idx<motif2_len; ++idx) {
        fake_dssp[idx+motif2_insert_pos] = 'X';
    }

    // ss prefix
    char pre_ss = dssp[0];
    Size pre_idx(0);
    std::vector<Size> span_len;
    std::vector<char> span_ss;
    for(Size idx=1; idx<length; ++idx) {
        if(dssp[idx] != pre_ss && pre_ss !=' ') {
            span_len.push_back(idx-pre_idx);
            span_ss.push_back(pre_ss);
            pre_idx = idx;
            pre_ss = dssp[idx];
        }
    }
    span_len.push_back(length-pre_idx);
    span_ss.push_back(pre_ss);
    std::stringstream prefix;
    for(Size idx=0; idx<span_len.size(); idx++) {
        prefix << span_ss[idx] << span_len[idx];
    }
    dssp_prefix = prefix.str();



    return pass;

}

bool random_motif_bundle_dssp_GFP(Size length,std::string motif_ss, Size num_helix, Size min_helix_len, Size max_helix_len, Size min_loop_len, Size max_loop_len, 
                                    Size GFP_upper_helix_max_len,
                                    Size GFP_upper_helix_min_len,
                                    Size GFP_lower_loop_max_len,
                                    Size GFP_lower_loop_min_len,
                                    Size GFP_lower_helix_max_len,
                                    Size GFP_lower_helix_min_len,
                                    std::string & dssp_prefix, std::string & dssp,Size & insert_pos)
{

    bool pass = true;
    std::vector<Size> loops(num_helix-1, 0);
    std::vector<Size> helices(num_helix-1, 0);
    Size max_tries = 100000;


    Size ihelix_fluorophore;
    Size fluorophore_helix_len;
    Size GFP_upper_helix_len;
    Size GFP_lower_loop_len;
    Size GFP_lower_helix_len;

    while(max_tries--) {
        Size res_len = length;

        ihelix_fluorophore = random_int(1, num_helix-2);
        GFP_upper_helix_len = random_int(GFP_upper_helix_min_len, GFP_upper_helix_max_len);
        GFP_lower_loop_len  = random_int(GFP_lower_loop_min_len , GFP_lower_loop_max_len);
        GFP_lower_helix_len = random_int(GFP_lower_helix_min_len, GFP_lower_helix_max_len);
        fluorophore_helix_len = GFP_upper_helix_len + motif_ss.length() + GFP_lower_loop_len + GFP_lower_helix_len;

        res_len -= fluorophore_helix_len;

        // loop
        for(Size idx=0; idx<loops.size(); ++idx) {
            Size this_loop_len = random_int(min_loop_len, max_loop_len);
            loops[idx] = this_loop_len;
            res_len -= this_loop_len;
        }

        // helix
        for(Size idx=0; idx<helices.size(); ++idx) helices[idx] = 0;
        while(res_len--) {
            Size idx = random_int(0, helices.size()-1);
            helices[idx] += 1;
        }

        pass = true;
        if (helices[0] < 16 || helices[num_helix-2] < 16) pass = false;
        for(Size idx=0; idx<helices.size(); ++idx) { if(helices[idx]<min_helix_len || helices[idx]>max_helix_len) pass = false; }
        if(pass) break;
    }

    helices.insert(helices.begin()+ihelix_fluorophore, fluorophore_helix_len);
    

    dssp_prefix = "";
    dssp = "";
    std::stringstream prefix;
    for(Size idx=0; idx<num_helix-1; idx++) {
        if(idx==ihelix_fluorophore) {
            dssp += std::string(GFP_upper_helix_len, 'H') + motif_ss + std::string(GFP_lower_loop_len, 'L') + std::string(GFP_lower_helix_len, 'H');
            prefix << "H" << GFP_upper_helix_len << "X" << motif_ss.length() << "L" << GFP_lower_loop_len << "H" << GFP_lower_helix_len;
        } else {
            dssp += std::string(helices[idx], 'H');
            prefix << "H" << helices[idx];
        }
        dssp += std::string(loops[idx], 'L');
        prefix << "L" << loops[idx];
    }
    dssp += std::string(helices[num_helix-1], 'H');
    prefix << "H" << helices[num_helix-1];

    dssp_prefix = prefix.str();

    insert_pos = 0;
    for(Size idx=0; idx<ihelix_fluorophore; ++idx) {
        insert_pos += helices[idx] + loops[idx]; 
    }
    insert_pos += GFP_upper_helix_len + 1;
    
    return pass;

}

bool random_motif_bundle_Nter_extension_dssp(Size length,std::string motif_ss,Size num_helix, Size min_helix_len, Size max_helix_len, Size min_loop_len, Size max_loop_len, std::string & dssp_prefix, std::string & dssp)
{

    bool pass = true;
    std::vector<Size> loops(num_helix, 0);
    std::vector<Size> helices(num_helix, 0);

    Size motif_len = motif_ss.length();

    Size max_tries = 100000;

    // Size min_terminal_helix_len = int( (length - (num_helix-1)*4)/num_helix );
    // if (min_terminal_helix_len>14) min_terminal_helix_len = 14;

    while(max_tries--) {
        Size res_len = length - motif_len;
        
        for(Size idx=0; idx<loops.size(); ++idx) {
            Size this_loop_len = random_int(min_loop_len, max_loop_len);
            loops[idx] = this_loop_len;
            res_len -= this_loop_len;
        }

        for(Size idx=0; idx<helices.size(); ++idx) helices[idx] = 0;
        while(res_len--) {
            Size idx = random_int(0, num_helix-1);
            helices[idx] += 1;
        }

        pass = true;
        // if (helices[0] < 14 || helices[num_helix-1] < 14) pass = false;
        if(helices[num_helix-1]<14) pass = false;
        for(Size idx=0; idx<helices.size(); ++idx) { if(helices[idx]<min_helix_len || helices[idx]>max_helix_len) pass = false; }
        if(pass) break;
    }


    dssp = motif_ss;
    std::stringstream prefix;
    prefix << "M" << motif_len;

    for(Size idx=0; idx<num_helix; idx++) {
        dssp += std::string(loops[idx], 'L') + std::string(helices[idx], 'H');
        prefix << "L" << loops[idx] << "H" << helices[idx];
    }

    dssp_prefix = prefix.str();

    return pass;

}


bool random_repeat_bundle_dssp(Size length, Size num_repeats, Size num_helix, Size min_helix_len, Size max_helix_len, Size min_loop_len, Size max_loop_len, std::string & dssp_prefix, std::string & dssp)
{

    assert( length % num_repeats == 0 );

    bool pass = true;
    std::vector<Size> loops(num_helix, 0);
    std::vector<Size> helices(num_helix, 0);

    Size max_tries = 100000;
    while(max_tries--) {
        Size res_len = length / num_repeats;
        
        for(Size idx=0; idx<loops.size(); ++idx) {
            Size this_loop_len = random_int(min_loop_len, max_loop_len);
            loops[idx] = this_loop_len;
            res_len -= this_loop_len;
        }

        for(Size idx=0; idx<helices.size(); ++idx) helices[idx] = 0;
        while(res_len--) {
            Size idx = random_int(0, num_helix-1);
            helices[idx] += 1;
        }

        pass = true;

	// the length of the first/last helix can not be too small
        if (helices[0] < 14 || helices[num_helix-1] < 14) pass = false;

        for(Size idx=0; idx<helices.size(); ++idx) { if(helices[idx]<min_helix_len || helices[idx]>max_helix_len) pass = false; }
        
        if(pass) break;
    }

    dssp_prefix = "";
    dssp = "";
    std::stringstream prefix;
    prefix<< "R" << num_repeats;
    for(Size idx=0; idx<num_repeats; idx++) {
        for(Size jdx=0; jdx<num_helix; jdx++) {
            dssp += std::string(helices[jdx], 'H') + std::string(loops[jdx], 'L');
        }
    }
    for(Size jdx=0; jdx<num_helix; jdx++) {
        prefix << "H" << helices[jdx] << "L" << loops[jdx];
    }

    dssp_prefix = prefix.str();

    return pass;
}

Real sidechain_neighbors(scene::Pose const & pose, bool including_intra_chain, bool including_inter_chain)
{
    auto per_res_scnb = per_res_sidechain_neighbors(pose, including_intra_chain, including_inter_chain);
    return per_res_scnb.sum() / pose.size();
}

Eigen::Array<Real, Eigen::Dynamic, Eigen::Dynamic>
per_res_sidechain_neighbors(scene::Pose const & pose, bool including_intra_chain, bool including_inter_chain)
{

    const Real dist_midpoint(9.0), dist_exponent(1.0), angle_shift_factor(0.5), angle_exponent(2.0);
    
    Size nres = pose.size();
    Size nchains = pose.num_chains();
    Size total_res = nres * nchains;

    Eigen::Matrix<Real, Eigen::Dynamic, 3> N(total_res, 3), 
                                           CA(total_res, 3),
                                           C(total_res, 3), 
                                           CB(total_res, 3), 
                                           N_CA(total_res, 3), 
                                           CA_C(total_res, 3), 
                                           CA_CB(total_res, 3);
    Eigen::Matrix<Real, Eigen::Dynamic, Eigen::Dynamic> dist_map(nres, total_res), angle_map(nres, total_res);

    Real total_score=0;

    for(Size ichain=1; ichain<=nchains; ++ichain) {
        for(Size idx=1; idx<=nres; ++idx) {
            Size global_idx = (ichain-1) * nres + idx - 1;
            CA.row(global_idx) = pose.xyz(idx, ATOM_CA, ichain);
            C.row(global_idx)  = pose.xyz(idx, ATOM_C, ichain);
            N.row(global_idx)  = pose.xyz(idx, ATOM_N, ichain);
        }
    }

    N_CA = CA - N;
    CA_C = C  - CA;
    for(Size idx=0; idx<total_res; ++idx)
        CA_CB.row(idx) = -0.58273431*N_CA.row(idx).cross(CA_C.row(idx)) + 0.56802827*N_CA.row(idx) - 0.54067466*CA_C.row(idx);
    CB = CA_CB + CA;
    CA_CB.rowwise().normalize();

    // intra_chain / inter_chain
    Size start_res(0), end_res(total_res);
    if (!including_intra_chain) start_res = nres;
    if(!including_inter_chain)  end_res   = nres;

    for(Size idx=0; idx<nres; ++idx){
        for(Size jdx=0; jdx<total_res; ++jdx){
            // necessory ????
            // not efficient, not smart
            angle_map(idx, jdx) = 0.0;
            dist_map(idx, jdx)  = 0.0;
        }
    }

    for(Size idx=0; idx<nres; ++idx)
    {
        for(Size jdx=start_res; jdx<end_res; ++jdx)
        {
            if(idx == jdx) continue;

            dist_map(idx, jdx) = 1.0/(1.0+std::exp(dist_exponent*((CB.row(jdx)-CB.row(idx)).norm() - dist_midpoint)));

            angle_map(idx, jdx) = (CA_CB.row(idx).dot((CB.row(jdx)-CB.row(idx)).normalized()) + angle_shift_factor) / (1+angle_shift_factor);
            if (angle_map(idx, jdx) > 0.0) {
                angle_map(idx, jdx) = std::pow( angle_map(idx, jdx), angle_exponent);
            } else {
                angle_map(idx, jdx) = 0.0;
            }
            /*
            // This actually works pretty well
            // To make the code consistent with the python script       
            angle_map(idx, jdx) = CA_CB.row(idx).dot((CB.row(jdx)-CB.row(idx)).normalized());
            if(angle_map(idx,jdx) > 0.0) {
                angle_map(idx,jdx) = std::pow( (angle_map(idx,jdx) + angle_shift_factor) / (1+angle_shift_factor), 2.0 );
            }
            else {
                angle_map(idx, jdx) = 0.0;
            }
            */

        }
    }

    return (angle_map.array() * dist_map.array()).rowwise().sum();
    // return (angle_map.array() * dist_map.array()).sum() / nres;
}

Real ligand_neighbors(scene::Pose const & pose, Eigen::Matrix<Real, Eigen::Dynamic, 3> const & ligand)
{
    const Real dist_midpoint(5.0), dist_exponent(1.0), angle_shift_factor(0.5), angle_exponent(1.0);
    
    Size nres = pose.size();
    Size nchains = pose.num_chains();
    Size total_res = nres * nchains;

    Size ligand_natoms = ligand.rows();

    std::string const sequence = pose.sequence();

    Eigen::Matrix<Real, Eigen::Dynamic, 3> N(total_res, 3), 
                                           CA(total_res, 3),
                                           C(total_res, 3), 
                                           CB(total_res, 3),
                                           /*
                                           N_CA(total_res, 3), 
                                           CA_C(total_res, 3), */
                                           CA_CB(total_res, 3);
                                        
    Eigen::Matrix<Real, Eigen::Dynamic, Eigen::Dynamic> dist_map(total_res, ligand_natoms), angle_map(total_res, ligand_natoms);

    Real total_score=0;

    for(Size ichain=1; ichain<=nchains; ++ichain) {
        for(Size idx=1; idx<=nres; ++idx) {
            Size global_idx = (ichain-1) * nres + idx - 1;
            CA.row(global_idx) = pose.xyz(idx, ATOM_CA, ichain);
            C.row(global_idx)  = pose.xyz(idx, ATOM_C, ichain);
            N.row(global_idx)  = pose.xyz(idx, ATOM_N, ichain);

            if( sequence.at(idx-1) == 'G' ) {
                Vec N_CA = CA.row(global_idx) - N.row(global_idx);
                Vec CA_C = C.row(global_idx) - CA.row(global_idx);
                CB.row(global_idx) = -0.58273431*N_CA.cross(CA_C) + 0.56802827*N_CA - 0.54067466*CA_C;
                CB.row(global_idx) += CA.row(global_idx);
            } else {
                CB.row(global_idx) = pose.xyz(idx, ATOM_CB, ichain);
            }
        }
    }

    /*
    N_CA = CA - N;
    CA_C = C  - CA;
    for(Size idx=0; idx<total_res; ++idx)
        CA_CB.row(idx) = -0.58273431*N_CA.row(idx).cross(CA_C.row(idx)) + 0.56802827*N_CA.row(idx) - 0.54067466*CA_C.row(idx);
    CB = CA_CB + CA;
    */
    CA_CB = CB - CA;
    CA_CB.rowwise().normalize();


    // initialize to 0
    for(Size idx=0; idx<total_res; ++idx){
        for(Size jdx=0; jdx<ligand_natoms; ++jdx){
            // necessory ????
            // not efficient, not smart
            angle_map(idx, jdx) = 0.0;
            dist_map(idx, jdx)  = 0.0;
        }
    }

    // fill in the numbers
    for(Size idx=0; idx<total_res; ++idx)
    {
        for(Size jdx=0; jdx<ligand_natoms; ++jdx)
        {

            dist_map(idx, jdx) = 1.0/(1.0+std::exp(dist_exponent*((ligand.row(jdx)-CB.row(idx)).norm() - dist_midpoint)));

            angle_map(idx, jdx) = (CA_CB.row(idx).dot((ligand.row(jdx)-CB.row(idx)).normalized()) + angle_shift_factor) / (1+angle_shift_factor);
            if (angle_map(idx, jdx) > 0.0) {
                angle_map(idx, jdx) = std::pow( angle_map(idx, jdx), angle_exponent);
            } else {
                angle_map(idx, jdx) = 0.0;
            }
            /*
            // This actually works pretty well
            // To make the code consistent with the python script       
            angle_map(idx, jdx) = CA_CB.row(idx).dot((CB.row(jdx)-CB.row(idx)).normalized());
            if(angle_map(idx,jdx) > 0.0) {
                angle_map(idx,jdx) = std::pow( (angle_map(idx,jdx) + angle_shift_factor) / (1+angle_shift_factor), 2.0 );
            }
            else {
                angle_map(idx, jdx) = 0.0;
            }
            */

        }
    }


    return (angle_map.array() * dist_map.array()).sum() / ligand_natoms;
}

Real
get_dihedral(Vec const & a, Vec const & b, Vec const & c, Vec const & d)
{
    Vec b0 = -1.0*(b-a);
    Vec b1 = (c - b);
    Vec b2 = d - c;

    b1.normalize();

    Vec v = b0 - b0.dot(b1) * b1;
    Vec w = b2 - b2.dot(b1) * b1;

    return std::atan2(b1.cross(v).dot(w), v.dot(w));
}

Real
get_angle(Vec const & a, Vec const & b, Vec const & c)
{
    Vec v = a - b;
    v.normalize();
    Vec w = c - b;
    w.normalize();

    return std::acos(v.dot(w));
}


EigenXform
xform_from_3points(Vec const & a, Vec const & b, Vec const & c)
{
    EigenXform X;
    X.translation() = b;

    Vec e1 = (c - b).normalized();
    Vec e3 = e1.cross(a-b).normalized();
    Vec e2 = e3.cross(e1).normalized();
    X.matrix().col(0) = e1;
    X.matrix().col(1) = e2;
    X.matrix().col(2) = e3;

    return X;

}

void print_xform(EigenXform const & x)
{
    std::cout << x(0,0) << " " << x(0,1) << " " << x(0,2) << " " << x(0,3) << std::endl
              << x(1,0) << " " << x(1,1) << " " << x(1,2) << " " << x(1,3) << std::endl
              << x(2,0) << " " << x(2,1) << " " << x(2,2) << " " << x(2,3) << std::endl;
}

void print_vec(Vec const & v)
{
    std::cout << v(0) << " " << v(1) << " " << v(2) << std::endl;
}


void get_repeat_parameters_from_stubs(scene::Pose & pose,Real & rise_out, Real & radius_out, Real & omega_out,Vec & axis_out,Vec & axis_center,bool debug )
{
    using namespace Eigen;
    assert(pose.num_repeats()>1);
    assert(pose.size()%pose.num_repeats() ==0);
	const Real FLOAT_PRECISION = 1e-6;

    Size residue_offset = pose.size()/pose.num_repeats();
	EigenXform sym_operator =pose.stub(1+residue_offset)*pose.stub(1).inverse() ;
	Eigen::AngleAxis<Real> angle_axis(sym_operator.rotation());

	// axis ang
	Vec axis = angle_axis.axis();
	Real ang = angle_axis.angle();

	Real translation_along_axis_of_rotation = axis.dot(sym_operator.translation());
	Vec cen(0.0,0.0,0.0);
	// these points lie on a circle and the plane is perpendicular with the axis 
    Vec c_A = pose.conformation().center_vec(1,residue_offset);
    Vec c_B = pose.conformation().center_vec(residue_offset+1,residue_offset*2);
    Vec p0 = c_A;
	Vec p1 = sym_operator * p0;
	Vec p2 = sym_operator * p1;
    if(debug){
        char buf[128];
        Size anum(1), rnum(1);
        std::ofstream out("/home/chentong/wefold/build/test_point.pdb");
        Vec xyz = p0;
        snprintf(buf,128,"%s%5i %4s %3s %c%4i    %8.3f%8.3f%8.3f%6.2f%6.2f %11s\n",
                    "HETATM",
                    anum++,
                    "BURR",
                    "BUR",
                    'B',
                    rnum++,
                    xyz[0],xyz[1],xyz[2],
                    1.0,
                    1.0,
                    "B"
                );

        out << buf;
        xyz = p1;
        snprintf(buf,128,"%s%5i %4s %3s %c%4i    %8.3f%8.3f%8.3f%6.2f%6.2f %11s\n",
                    "HETATM",
                    anum++,
                    "BURR",
                    "BUR",
                    'B',
                    rnum++,
                    xyz[0],xyz[1],xyz[2],
                    1.0,
                    1.0,
                    "B"
                );

        out << buf;
        xyz = p2;
        snprintf(buf,128,"%s%5i %4s %3s %c%4i    %8.3f%8.3f%8.3f%6.2f%6.2f %11s\n",
                    "HETATM",
                    anum++,
                    "BURR",
                    "BUR",
                    'B',
                    rnum++,
                    xyz[0],xyz[1],xyz[2],
                    1.0,
                    1.0,
                    "B"
                );

        out << buf;
    }
	p1 -= axis * (p1-p0).dot(axis);
	p2 -= axis * (p2-p0).dot(axis);
	Real d = (p1-p0).norm();

    Real l = d / (2.0*std::tan(ang/2.0));
    Vec tocen = (p1-p0).normalized().cross(axis)*l;
    if(tocen.dot(p2-p1)<0.0){
        tocen = -tocen;
    }
    cen = (p0+p1)/2.0 + tocen;
    Real radius = (p0-cen).norm();
    radius_out =radius;
    omega_out = ang;
    rise_out = translation_along_axis_of_rotation;
    axis_out = axis;
    
    axis_center = cen + translation_along_axis_of_rotation * axis*(float(pose.num_repeats())/2.0 -0.5);
    return;
}
void get_repeat_parameters_from_coords(scene::Pose & pose, Real & rise_out, Real & radius_out, Real & omega_out,Vec & axis_out,Vec & axis_center)  {
	using Eigen::MatrixXd;
	using namespace Eigen;
	using namespace std;
	Size seg_size = pose.size()/pose.num_repeats();
    Size startAtRepeat = pose.num_repeats()/2;
	Size startResOffset = seg_size*(startAtRepeat-1);
	MatrixXf A(3, seg_size);
	MatrixXf B(3, seg_size);
	for ( Size i = 1; i <= seg_size; i++ ) {
		Size offset =i+startResOffset;
		Vec coord = pose.xyz(offset,ATOM_CA);
		A(0,i-1) = coord[0];
		A(1,i-1) = coord[1];
		A(2,i-1) = coord[2];

	}
	for ( Size i = seg_size+1; i <= seg_size*2 ; i++ ) {
		Size offset =i+startResOffset;
		Vec coord = pose.xyz(offset,ATOM_CA);
		B(0,i-1-seg_size) = coord[0];
		B(1,i-1-seg_size) = coord[1];
		B(2,i-1-seg_size) = coord[2];
	}
	Vector3f c_A(A.row(0).mean(), A.row(1).mean(), A.row(2).mean());
	Vector3f c_B(B.row(0).mean(), B.row(1).mean(), B.row(2).mean());
	MatrixXf x_A(A.rows(),A.cols());
	MatrixXf x_B(B.rows(),B.cols());
	for ( int i=0; i<A.cols(); i++ ) {
		x_A.col(i)=A.col(i)-c_A;
		x_B.col(i)=B.col(i)-c_B;
	}
	Matrix3f cov= (x_B * x_A.transpose()) / x_A.cols();
	JacobiSVD<MatrixXf> svd(cov, ComputeFullU | ComputeFullV);

	Matrix3f Rt=svd.matrixU() * svd.matrixV().transpose();
	Matrix3f R;
	R<< 1,0,0, 0,1,0, 0,0,Rt.determinant();
	Matrix3f H= svd.matrixU() * R * svd.matrixV().transpose();

	Real acos_fix_tmp = (H.trace()-1)/2;
	if ( acos_fix_tmp <= -1.0 ) {
		acos_fix_tmp = -1.0;
	}
	if ( acos_fix_tmp >= 1.0 ) {
		acos_fix_tmp = 1.0;
	}
	Real omega = acos(acos_fix_tmp);
	Matrix3f I = Matrix3f::Identity();
	Matrix3f N = 0.5*(H+H.transpose()) - cos(omega)*I;


	Vector3f hN(0,0,0);//initializiation for the testing server
	Real scalar = 0;
	Real max_scalar = -10000000;
	for ( Size i = 0; i<=2; i++ ) {
		scalar = N.col(i).norm();
		if ( scalar > max_scalar ) {
			max_scalar = scalar;
			hN = N.col(i)/N.col(i).norm();
		}
	}

	Real sin_omega = (H(1,0)-H(0,1)) / (2*hN(2));
	if ( sin_omega < 0 ) hN = -1 * hN;

	Vector3f t = c_B - H*c_A;
	Real L = t.dot(hN) ;
	Real rise=abs(L);

	Matrix3f Ncross;
	Ncross << 0,-1*hN(2),hN(1), hN(2),0,-1*hN(0), -1*hN(1),hN(0),0 ;
	Matrix3f R0t= (1-cos(omega))*I - sin(omega)*Ncross;
	Vector3f R0 = R0t.inverse() * (t-L*hN);
	Vector3f pA= (c_A-R0)-(hN*(hN.dot(c_A-R0)));
	Vector3f pB= (c_B-R0)-(hN*(hN.dot(c_B-R0)));


	Real direction = L * hN.dot(pA.cross(pB));
	Real radius=pA.norm();
	rise_out = rise;
	radius_out = radius;
	omega_out = omega;
    axis_out = hN;
    axis_center = R0+(hN*(hN.dot(c_A-R0))) + hN*(float(pose.num_repeats())/2+0.5-(pose.num_repeats()/2));
}

Real rmsd_no_super(const scene::Pose & pose1, const scene::Pose & pose2, Size start_res, Size end_res, bool CA_only)
{

    if( end_res == -1 ) end_res = pose1.size();
    assert( pose1.size() == pose2.size() && start_res < end_res && end_res <= pose1.size() );

    std::vector<ATOM_TYPE> atom_types = {ATOM_CA, ATOM_N, ATOM_C};

    Real dev(0.0);
    Size counts(0);
    for(Size ires=start_res; ires<=end_res; ++ires) {
        for( Size iatom=0; iatom<3; ++iatom ) {

            Vec  iatom_xyz = pose1.xyz(ires, atom_types[iatom]);
            Vec  jatom_xyz = pose2.xyz(ires, atom_types[iatom]);
            dev += (iatom_xyz - jatom_xyz).squaredNorm();
            ++counts;

            if(CA_only) break;
        }
    }

    return std::sqrt(dev / counts);
}

std::string get_dssp_from_pose(scene::Pose pose ,Size len,bool reduece_ss){
    std::string dssp;
    std::vector<std::vector<std::vector<double>>> coords;
    std::vector<ATOM_TYPE> atom_types{ATOM_N,ATOM_CA,ATOM_C,ATOM_O};
    for(size_t ires=1;ires<=(len==-1?pose.size():len);ires++){
        coords.push_back(std::vector<std::vector<double>>());
        for(size_t j=0;j<4;j++){
            Vec atom_vec = pose.xyz(ires,atom_types[j]);
            coords[ires-1].push_back(std::vector<double>{atom_vec.x(),atom_vec.y(),atom_vec.z()});
        }
    }
    MProtein a;
    a.ReadCoords(coords);
    a.CalculateSecondaryStructure();
    std::vector<const MResidue*> residues;

    for (const MChain* chain:a.GetChains())
    {
        for (const MResidue* residue:chain->GetResidues())
        residues.push_back(residue);
    }

    for (const MResidue* res: residues){
        char ss;
        const MResidue residue = *res;
        if(reduece_ss){
            switch (residue.GetSecondaryStructure())
            {
                case alphahelix:  ss = 'H'; break;
                case betabridge:  ss = 'E'; break;
                case strand:    ss = 'E'; break;
                case helix_3:    ss = 'H'; break;
                case helix_5:    ss = 'H'; break;
                case turn:      ss = 'L'; break;
                case bend:      ss = 'L'; break;
                case loop:      ss = 'L'; break;
            }
        }else{
            switch (residue.GetSecondaryStructure())
            {
                case alphahelix:  ss = 'H'; break;
                case betabridge:  ss = 'B'; break;
                case strand:    ss = 'E'; break;
                case helix_3:    ss = 'G'; break;
                case helix_5:    ss = 'I'; break;
                case turn:      ss = 'T'; break;
                case bend:      ss = 'S'; break;
                case loop:      ss = ' '; break;
            }
        }
        std::string NHO[2], ONH[2];
        int64 nNHO[2], nONH[2];
        const HBond* acceptors = residue.Acceptor();
        const HBond* donors = residue.Donor();
        for (uint32 i = 0; i < 2; ++i)
        {
            NHO[i] = ONH[i] = "0, 0.0";
            nNHO[i] = nONH[i] = 0;

            if (acceptors[i].residue != nullptr)
            {
            nNHO[i] = acceptors[i].residue->GetNumber() - residue.GetNumber();
            NHO[i] =acceptors[i].energy;
            }

            if (donors[i].residue != nullptr)
            {
            nONH[i] = donors[i].residue->GetNumber() - residue.GetNumber();
            ONH[i] = donors[i].energy;
            }
        }
        dssp.push_back(ss);
    }

    return dssp;
}

}
