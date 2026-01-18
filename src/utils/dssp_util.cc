#include "utils/dssp_util.hh"
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

std::string normalize_bundle_dssp(std::string dssp_str) {
    std::string ss = dssp_str;
    Size ss_len = ss.length();
    if( ss[2] == 'H' && ss[3] == 'H' && ss[4] == 'H' && ss[5] == 'H' && 
        (ss[0] != 'H' || ss[1] != 'H' ) ) {
        ss[0] = 'H';
        ss[1] = 'H';
    }
    if( ss[ss_len-3] == 'H' && ss[ss_len-4] == 'H' && ss[ss_len-5] == 'H' && ss[ss_len-6] == 'H' && 
        (ss[ss_len-2] != 'H' || ss[ss_len-1] != 'H' ) ) {
        ss[ss_len-1] = 'H';
        ss[ss_len-2] = 'H';
    }
    return ss;
}



}
