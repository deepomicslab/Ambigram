#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>
#include <algorithm>
#include <getopt.h>
#include <map>

#include <htslib/sam.h>
#include <htslib/tbx.h>
//deprecated
using namespace std;

const map<string, int> SV_TYPE = {
        {"dup",   0},
        {"del",   1},
        {"trans", 2},
        {"inv",   3},
        {"ins",   4}
};
typedef struct Segment {
    int id;
    string chrom;
    int start;
    int end;
    double depth;
    bool is_inv = false;
    bool is_ins = false;
} seg_t;

typedef struct SequenceMap {
    int s_start, s_end;
    int r_start, r_end;
    char inv_indicator;
} seq_map_t;

typedef struct Junction {
    seg_t *left;
    seg_t *right;
    int support;
} junc_t;

/*
 * @sv_type:
 *      dup: 0
 *      del: 1
 *      trans: 2
 *      inv: 3
 *      ins: 4
*/
vector<seq_map_t *> read_data(char const *fn, int sv_type = -1) {
    vector<seq_map_t *> seq_maps;
    ifstream inf(fn);
    string line;
    stringstream ss;
    int s_start, s_end;
    int r_start, r_end;
    char inv_indicator;
    bool is_inv = false;
    bool is_ins = false;
    getline(inf, line); // skip header
    int ins_id = 1;
    while (getline(inf, line)) {
        ss << line;
        seq_map_t *row;
        switch (sv_type) {
            case 0:
            case 1:
            case 2:
                ss >> s_start >> s_end >> r_start >> r_end;
                break;
            case 3:
                ss >> s_start >> s_end >> r_start >> r_end >> inv_indicator;
                break;
            case 4:
                break;
        }
        row = new seq_map_t{s_start, s_end, r_start, r_end, inv_indicator};
        seq_maps.push_back(row);
        ss.clear();
    }
    inf.close();
    return seq_maps;
}

vector<seg_t *> read_segs(char const *fn) {
    vector<seg_t *> segs;

    ifstream inf(fn);
    string line;
    stringstream ss;
    int id, start, end;
    double count;
    string chrom;
    getline(inf, line);
    while (getline(inf, line)) {
        ss << line;
        ss >> id >> chrom >> start >> end >> count;
        seg_t *seg = new seg_t{id, chrom, start, end, 0};
        segs.push_back(seg);
        ss.clear();
    }
    inf.close();
    return segs;
}


vector<seg_t *> get_seg_seq(vector<seq_map_t *> &seq_maps, vector<seg_t *> &segs) {
    vector<seg_t *> seg_seq;
    vector<seg_t *>::iterator itr;
    for (seq_map_t *row : seq_maps) {
        itr = find_if(segs.begin(), segs.end(), [row](seg_t *s) { return s->start == row->r_start; });
        (*itr)->is_inv = row->inv_indicator == 'I';
        seg_seq.push_back(*itr);
    }
    return seg_seq;
}

vector<junc_t *> get_juncs(vector<seg_t *> &seg_seq) {
    vector<junc_t *> juncs;
    seg_t *left;
    seg_t *right;
    bool found;
    for (int i = 0; i < seg_seq.size() - 1; i++) {
        left = seg_seq[i];
        right = seg_seq[i + 1];
        found = false;
        for (junc_t *j : juncs) {
            if (j->left == left && j->right == right) {
                found = true;
                break;
            }
        }
        if (!found) {
            junc_t *junc = new junc_t{left, right, 0};
            juncs.push_back(junc);
        }
    }
    return juncs;
}

int get_overlap(bam1_t *aln, int start, int end) {
    uint32_t *cigar = bam_get_cigar(aln);
    int aln_start = aln->core.pos;
    int aln_end = aln_start + bam_cigar2rlen(aln->core.n_cigar, cigar) - 1;
    int overlap_start, overlap_end;
    if (aln_start > start) {
        overlap_start = aln_start;
    } else {
        overlap_start = start;
    }
    if (aln_end > end) {
        overlap_end = end;
    } else {
        overlap_end = aln_end;
    }
    // cout << "overlap start: " << overlap_start << " overlap end: " << overlap_end << endl;
    return overlap_end - overlap_start + 1;
}

void count_support(char const *bam_fn, vector<junc_t *> &juncs) {
    samFile *bam1 = sam_open(bam_fn, "rb");
    samFile *bam2 = sam_open(bam_fn, "rb");
    bam_hdr_t *header1 = sam_hdr_read(bam1);
    bam_hdr_t *header2 = sam_hdr_read(bam2);
    hts_idx_t *idx1 = sam_index_load(bam1, bam_fn);
    hts_idx_t *idx2 = sam_index_load(bam2, bam_fn);
    bam1_t *aln1 = bam_init1();
    bam1_t *aln2 = bam_init1();
    int cid1 = sam_hdr_name2tid(header1, "BOR");
    int cid2;
    hts_itr_t *iter1 = NULL;
    hts_itr_t *iter2 = NULL;
    uint8_t *aux;
    char *aux_token;
    int aux_pos;
    int left_pos, right_pos;

    for (junc_t *junc : juncs) {
        if (junc->left->is_inv) {
            left_pos = junc->left->start;
            iter1 = sam_itr_queryi(idx1, cid1, left_pos - 1, left_pos);
        } else {
            left_pos = junc->left->end;
            iter1 = sam_itr_queryi(idx1, cid1, left_pos - 1 - 1, left_pos - 1);
        }
        if (junc->right->is_inv) {
            right_pos = junc->right->end;
        } else {
            right_pos = junc->right->start;
        }
        // cout << "left: " << junc->left->id << " " << left_pos << ", " 
        //      << "right: " << junc->right->id << " " << right_pos << endl;
        // iter1 = sam_itr_queryi(idx1, cid1, left_pos - 1 - 1, left_pos - 1);
        if (left_pos == right_pos) {
            while (sam_itr_next(bam1, iter1, aln1) >= 0) {
                if (get_overlap(aln1, left_pos - 1 - 9 - 1, right_pos - 1 + 10) >= 20) {
                    junc->support++;
                }
            }
        } else {
            while (sam_itr_next(bam1, iter1, aln1) >= 0) {
                aux = bam_aux_get(aln1, "SA");
                if (aux != NULL) {
                    char *sa_chrom = strtok(bam_aux2Z(aux), ",");
                    aux_token = strtok(NULL, ",");
                    aux_pos = atoi(aux_token);
                    cid2 = sam_hdr_name2tid(header2, sa_chrom);
                    iter2 = sam_itr_queryi(idx2, cid2, aux_pos - 1, aux_pos);
                    while (sam_itr_next(bam2, iter2, aln2) >= 0) {
                        if (strcmp(bam_get_qname(aln2), bam_get_qname(aln1)) == 0) {
                            if (bam_aux_get(aln2, "SA") == NULL) continue;
                            int left_overlap, right_overlap;
                            if (junc->left->is_inv) {
                                left_overlap = get_overlap(aln1, left_pos - 1, left_pos - 1 + 10);
                            } else {
                                left_overlap = get_overlap(aln1, left_pos - 9 - 1 - 1, left_pos - 1);
                            }
                            if (junc->right->is_inv) {
                                right_overlap = get_overlap(aln2, right_pos - 9 - 1 - 1, right_pos - 1);
                                // cout << right_pos - 9 - 1 - 1 << " " << right_pos - 1 << endl;
                            } else {
                                right_overlap = get_overlap(aln2, right_pos - 1, right_pos - 1 + 10);
                            }
                            // if (left_overlap < 0 || right_overlap < 0) {
                            //     cout << "left: " << junc->left->id << " " << left_pos << ", " 
                            //          << "right: " << junc->right->id << " " << right_pos << endl;
                            //     cout << bam_get_qname(aln2) << endl;
                            // cout << left_overlap << " " << right_overlap << endl;
                            // }
                            if (left_overlap >= 10 && right_overlap >= 10) {
                                junc->support++;
                            }
                        }
                    }
                }
            }
        }
    }

    bam_destroy1(aln1);
    bam_destroy1(aln2);
    hts_idx_destroy(idx1);
    hts_idx_destroy(idx2);
    bam_hdr_destroy(header1);
    bam_hdr_destroy(header2);
    sam_close(bam1);
    sam_close(bam2);
}

int get_avg_depth(char const *depth_fn, string chrom, int start, int end) {
    htsFile *depth = hts_open(depth_fn, "rb");
    tbx_t *idx = tbx_index_load(depth_fn);
    int cid = tbx_name2id(idx, chrom.c_str());
    hts_itr_t *iter = tbx_itr_queryi(idx, cid, start - 1, end);
    kstring_t str = {0, 0, NULL};

    char *token;
    int p, d;
    int tot_depth = 0;
    while (tbx_itr_next(depth, idx, iter, &str) >= 0) {
        token = strtok(str.s, "\t ");
        p = atoi(strtok(NULL, "\t "));
        d = atoi(strtok(NULL, "\t "));
        tot_depth += d;
    }

    hts_itr_destroy(iter);
    tbx_destroy(idx);
    hts_close(depth);
    return tot_depth * 1.0 / (end - start + 1);
}

void get_avg_depth(vector<seg_t *> &segs, char const *depth_fn) {
    htsFile *depth = hts_open(depth_fn, "rb");
    tbx_t *idx = tbx_index_load(depth_fn);
    int cid;
    hts_itr_t *iter = NULL;
    for (seg_t *seg: segs) {
        cid = tbx_name2id(idx, seg->chrom.c_str());
        iter = tbx_itr_queryi(idx, cid, seg->start - 1, seg->end - 1);
        char *token;
        int p, d;
        int tot_depth = 0;
        kstring_t str = {0, 0, NULL};
        while (tbx_itr_next(depth, idx, iter, &str) >= 0) {
            token = strtok(str.s, "\t ");
            p = atoi(strtok(NULL, "\t "));
            d = atoi(strtok(NULL, "\t "));
            tot_depth += d;
        }
        seg->depth = tot_depth * 1.0 / (seg->end - seg->start + 1);
    }
    hts_itr_destroy(iter);
    tbx_destroy(idx);
    hts_close(depth);
}

int cal_avg_seg_depth(vector<seg_t *> &segs) {
    vector<int> depths;
    for (seg_t *seg: segs) {
        depths.push_back(seg->depth);
    }
    sort(depths.begin(), depths.end());
    return depths[depths.size() / 2];
}

int cal_avg_junc_depth(vector<junc_t *> &juncs) {
    vector<int> depths;
    for (junc_t *junc: juncs) {
        depths.push_back(junc->support);
    }
    sort(depths.begin(), depths.end());
    return depths[depths.size() / 2];
}

void generate_lh(char const *fn, vector<seg_t *> &segs, vector<junc_t *> &juncs, char const *samplename) {
    ofstream outf(fn);

    outf << "SAMPLE " << samplename << endl
         << "AVG_SEG_DP " << cal_avg_seg_depth(segs) << endl
         << "AVG_JUNC_DP " << cal_avg_junc_depth(juncs) << endl
         << "PURITY 1" << endl
         << "AVG_PLOIDY 1" << endl
         << "PLOIDY 1" << endl
         << "SOURCE H:1" << endl
         << "SINK H:" << segs.back()->id << endl;

    for (seg_t *seg : segs) {
        outf << "SEG H:" << seg->id << ":" << seg->chrom << ":"
             << seg->start << ":" << seg->end << " " << seg->depth << " " << -1 << endl;
    }
    for (junc_t *junc : juncs) {
        outf << "JUNC H:" << junc->left->id << ":+ H:"
             << junc->right->id << ":+ " << junc->support << " -1 U B" << endl;
    }

    outf.close();
}

void write_juncs(char const *fn, vector<junc_t *> &juncs) {
    ofstream outf(fn);
    outf << "id_5p\tid_3p\tsupport" << endl;
    for (junc_t *junc : juncs) {
        outf << junc->left->id << "\t"
             << junc->right->id << "\t"
             << junc->support << endl;
    }
    outf.close();
}

void write_segs(char const *fn, vector<seg_t *> &segs) {
    ofstream outf(fn);
    outf << "id\tchrom\tstart\tend\tdepth\tis_inv\tis_ins" << endl;
    for (seg_t *seg : segs) {
        outf << seg->id << "\t"
             << seg->chrom << "\t"
             << seg->start << "\t"
             << seg->end << "\t"
             << seg->depth << "\t"
             << (seg->is_inv ? "True" : "False") << "\t"
             << (seg->is_ins ? "True" : "False") << endl;
    }
    outf.close();
}

int main(int argc, char *argv[]) {

    // * Input    
    char const *seg_fn, *data_fn;
    char const *bam_fn, *depth_fn;
    char const *sv_type;
    // * Output
    char const *seg_out_fn, *junc_out_fn;

    static struct option long_options[] = {
            {"seg_file",  required_argument, 0, 's'},
            {"data_file", required_argument, 0, 'j'},
            {"bam",       required_argument, 0, 'b'},
            {"depth",     required_argument, 0, 'd'},
            {"svtype",    required_argument, 0, 't'},

            {"seg_out",   required_argument, 0, 'S'},
            {"junc_out",  required_argument, 0, 'J'},
            {0, 0,                           0, 0}
    };

    int opt_idx = 0;
    int c = getopt_long(argc, argv, "s:j:b:d:t:S:J:", long_options, &opt_idx);
    while (c > 0) {
        switch (c) {
            case 's':
                seg_fn = optarg;
                cout << "test" << endl;
                break;
            case 'j':
                data_fn = optarg;
                break;
            case 'b':
                bam_fn = optarg;
                break;
            case 'd':
                depth_fn = optarg;
                break;
            case 't':
                sv_type = optarg;
                break;
            case 'S':
                seg_out_fn = optarg;
                break;
            case 'J':
                junc_out_fn = optarg;
                break;
            default:
                abort();
        }
        c = getopt_long(argc, argv, "s:j:b:d:t:S:J:", long_options, &opt_idx);
    }
    cout << SV_TYPE.size() << endl;
    cout << "seg_file: " << seg_fn << endl
         << "data_file: " << data_fn << endl
         << "bam: " << bam_fn << endl
         << "depth: " << depth_fn << endl
         << "svtype: " << sv_type << " " << SV_TYPE.at(sv_type) << endl
         << "seg_out: " << seg_out_fn << endl
         << "junc_out_fn: " << junc_out_fn << endl;

    vector<seq_map_t *> seq_maps = read_data(data_fn, SV_TYPE.at(sv_type));
    vector<seg_t *> segs = read_segs(seg_fn);
    vector<seg_t *> seg_seq = get_seg_seq(seq_maps, segs);
    vector<junc_t *> juncs = get_juncs(seg_seq);

    get_avg_depth(segs, depth_fn);
    count_support(bam_fn, juncs);
    write_segs(seg_out_fn, segs);
    write_juncs(junc_out_fn, juncs);
    return 0;
}
