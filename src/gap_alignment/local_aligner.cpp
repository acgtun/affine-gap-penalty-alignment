#include "option.hpp"
#include "evalue.hpp"
#include "bio_util.hpp"
#include "fasta_file.hpp"
#include "local_alignment.hpp"

#include <vector>
#include <string>
#include <fstream>

using std::vector;
using std::string;
using std::ofstream;

using local_alignment::LocalAlignment;

/* Output the alingment results for one query to fout */
void DisplayResults(const string& query_name, const string& database_file,
                    const vector<M8Results>& aligned_results, const int& outfmt,
                    ofstream& fout) {
  if (outfmt == 7) {
    fout << "# S3 1.0.0 Jan, 2015" << endl;
    fout << "# Query: " << query_name << endl;
    fout << "# Database: " << database_file << endl;
    fout
        << "# Fields: query id, subject id, % identity, alignment length, mismatches, "
            "gap opens, q. start, q. end, s. start, s. end, evalue, bit score"
        << endl;
    fout << "# " << aligned_results.size() << " hits found" << endl;
  }
  for (uint32_t i = 0; i < aligned_results.size(); i++) {
    fout << query_name << "\t" << aligned_results[i].protein_name << "\t"
        << aligned_results[i].identity << "\t" << aligned_results[i].aligned_len
        << "\t" << aligned_results[i].mismatch << "\t"
        << aligned_results[i].gap_open << "\t" << aligned_results[i].qs << "\t"
        << aligned_results[i].qe << "\t" << aligned_results[i].ps << "\t"
        << aligned_results[i].pe << "\t" << aligned_results[i].evalue << "\t"
        << aligned_results[i].bit_score << endl;
  }
}

int main(int argc, const char* argv[]) {
  InitProgram(argc, argv);
  string database_file, query_file, output_file;
  int outfmt;
  Option::GetOption("-d", database_file);
  Option::GetOption("-q", query_file);
  Option::GetOption("-f", outfmt, 7);
  Option::GetOption("-o", output_file);

  FastaFile database(database_file);
  FastaFile queries(query_file);

  Evalue evalue(database.num_of_characters, database.num_of_sequences);
  LocalAlignment local_align(&evalue, queries.max_sequence_length,
                             database.max_sequence_length);

  ofstream fout(output_file.c_str());
  for (uint32_t i = 0; i < queries.num_of_sequences; ++i) {
    evalue.UpdateValues(strlen(queries.sequences[i]));
    vector<M8Results> aligned_results;
    M8Results res;
    for (uint32_t j = 0; j < database.num_of_sequences; ++j) {
      if (local_align.RunLocalAlignment(queries.sequences[i],
                                        database.sequences[j],
                                        strlen(queries.sequences[i]),
                                        strlen(database.sequences[j]), res)) {
        res.protein_name = database.sequences_names[j];
        aligned_results.push_back(res);
      }
    }
    DisplayResults(queries.sequences_names[i], database_file,
                   aligned_results, outfmt, fout);
  }

  fout.close();

  return 0;
}
