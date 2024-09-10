#include <RtypesCore.h>
#include <TBrowser.h>
#include <TCanvas.h>
#include <TF1.h>
#include <TFile.h>
#include <TFitResult.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TPad.h>
#include <TPaveStats.h>
#include <TProfile.h>
#include <TROOT.h>
#include <TString.h>
#include <TStyle.h>
#include <TTree.h>
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdio>
#include <ctime>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <stdint.h>
#include <stdlib.h>
#include <string>
#include <sys/stat.h>
#include <unistd.h>
#include <vector>

TString filename_data;
std::vector<std::string> file_names;
int vcnt[128][32];
// averages channel in the neighborhood for dynamic threshold
std::vector<float> neigh_even;
std::vector<float> neigh_odd;

void analyze() {
  //! Variables for reading, analyzing data and displaying histograms
  int ch_min = 0;
  int ch_max = 128;
  int ch_step = 1;

  int d_min = 0;
  int d_max = 31;
  int d_step = 1;

  filename_data = (filename_data + ".txt");
  std::ifstream scanfile;
  scanfile.open(filename_data);
  int fch;
  std::string line;
  std::string ele;

  int count = 0;
  int ch_cnt = 0;
  while (!scanfile.eof()) {
    std::getline(scanfile, line);
    std::istringstream iss(line);
    if (!line.empty()) {
      iss >> ele;
      iss >> ele;
      count = 0;
      for (int d = d_min; d < d_max + 1; d++) {
        iss >> count;
        vcnt[ch_cnt][d] = count;
      }
      ch_cnt++;
    }
  }

  scanfile.close();
}

bool read_file_tests() {
  // ----------------------------------------------------------------------------------------------------------------
  TString file_list = "a";
  // ----------------------------------------------------------------------------------------------------------------
  std::ifstream myfile;
  myfile.open(file_list + ".txt");
  std::string line;
  std::string delimiter = ".txt";
  std::size_t pos = 0;
  std::string token;

  // reading list file and storing it in a vector
  while (!myfile.eof()) {
    std::getline(myfile, line);
    if (line.empty())
      continue;
    else {
      while ((pos = line.find(delimiter)) != std::string::npos) {
        token = line.substr(0, pos);
        // cout << token << endl;
        line.erase(0, pos + delimiter.length());
      }
      file_names.push_back(token);
    }
  }
  if (file_names.size() != 0)
    return false;
  return true;
}

// get_file_name from the file list
std::string get_file_name(int i) { return file_names[i]; }

// write the summary of broken channels for each file
void write_results(TString filename_str, std::vector<int> brk_ch) {
  std::ofstream myfile;
  myfile.open("conn_check_summary.txt", std::ios::out | std::ios::app);
  myfile << filename_str << ": ";
  if (brk_ch.size() != 0) {
    for (int i = 0; i < brk_ch.size(); i++) {
      myfile << brk_ch[i] << ",";
    }
    myfile << "\n";
  } else {
    myfile << "0"
           << "\n";
  }
  myfile.close();
}

/*  Check if a file exists
    returns  true if the file exists, else false */
bool file_exists(const std::string &filename) {
  struct stat buf;
  if (stat(filename.c_str(), &buf) != -1) {
    return true;
  }
  return false;
}
// write the parameters used to find the total number of broken channels
void write_parameters(float z_alpha, int n_neigh_ch,
                      std::size_t n_broken_channels) {
  /* z_alpha - percentage of neighbors' median for the threshold
   n_neigh_ch - number of channels in the neighborhood for dynamic threshold */
  std::ofstream my_file;
  if (!file_exists("conn_check_parameters.txt")) {
    my_file.open("conn_check_parameters.txt", std::ios::out);

    my_file << std::fixed << std::setw(20) << "z_alpha" << '\t' << std::fixed
            << std::setw(20) << "n_neigh_ch" << '\t' << std::fixed
            << std::setw(20) << "        " << '\t' << std::fixed
            << std::setw(20) << "        " << '\t' << std::fixed
            << std::setw(20) << "n_broken_channels" << '\n';

    my_file << std::fixed << std::setw(20) << z_alpha << '\t' << std::fixed
            << std::setw(20) << n_neigh_ch << '\t' << std::fixed
            << std::setw(20) << "        " << '\t' << std::fixed
            << std::setw(20) << "        " << '\t' << std::fixed
            << std::setw(20) << n_broken_channels << '\n';
  } else {
    my_file.open("conn_check_parameters.txt", std::ios::app);

    my_file << std::fixed << std::setw(20) << z_alpha << '\t' << std::fixed
            << std::setw(20) << n_neigh_ch << '\t' << std::fixed
            << std::setw(20) << "        " << '\t' << std::fixed
            << std::setw(20) << "        " << '\t' << std::fixed
            << std::setw(20) << n_broken_channels << '\n';
  }

  my_file.close();
}

template <typename T> T calculate_median(std::vector<T> v, std::size_t n) {
  // Sort the vector
  std::sort(v.begin(), v.end());
  // Check if the number of elements is odd
  if (n % 2 != 0) {
    return (double)v[n / 2];
  }
  // If the number of elements is even, return the average
  // of the two middle elements
  return (double)(v[(n - 1) / 2] + v[n / 2]) / 2.0;
}

// Write a vector of channels to the histogram
void fill_hist_with_ch(TH1F *hist, std::vector<int> channels, int *ch_hits) {
  if (channels.size() != 0) {
    for (int j = 0; j < channels.size(); j++) {
      hist->Fill(channels[j], ch_hits[channels[j]]);
    }
  }
}

// ------------------------ Main function ---------------------
int analysis_conn_check() {

  for (int ch = 0; ch < 128; ch++) {
    for (int d = 0; d < 32; d++) {
      vcnt[ch][31 - d] = 0;
    }
  }

  TString dir = "";
  TString filename_root = "";
  std::vector<int> broken_channels;

  // counter of broken channels accross all the files
  std::size_t n_broken_channels = 0;
  // percentage of neighbors' median for the threshold
  float z_alpha = 0.55;
  // number of channels in the neighborhood for dynamic threshold
  int n_neigh_ch = 25;

  // reading file with list of measurements files.
  read_file_tests();
  // analyze function is called for every file in the list. At the end,
  // everything is closed and root files can be accessed via TBrowser
  for (int i = 0; i < int(file_names.size()); i++) {
    filename_data = dir + get_file_name(i).c_str();
    filename_root = dir + get_file_name(i).c_str();
    // cout << "Filename_selected: " << filename_data << endl;
    analyze();
    // erase the previous channels
    broken_channels.clear();

    // Creating files
    TFile *file1 = new TFile(filename_root + ".root", "recreate");
    file1->cd();
    // Creating histograms
    TH1F h_ave_hits("h_hits", "h_hits", 128, 0, 128);
    TH1F h_broken_ch("h_broken_ch", "h_broken_ch", 128, 0, 128);

    TH1F h1_ave_even("h1_ave_even", "h1_ave_even", 100, 0, 400000);
    TH1F h1_ave_odd("h1_ave_odd", "h1_ave_odd", 100, 0,
                    400000); // access these channels
    TH1F h_ave_even("h_ave_even", "h_ave_even", 64, 0, 128);
    TH1F h_ave_odd("h_ave_odd", "h_ave_odd", 64, 1, 129);
    TH2F h2hits_amp("h2hits_amp", "h2hits_amp", 128, 0, 128, 31, 1, 31);

    // histograms storing thresholds for each channel, even and odd separately
    TH1F h_thr_even("h_thr_even", "h_thr_even", 128, 0, 128);
    TH1F h_thr_odd("h_thr_odd", "h_thr_odd", 128, 0, 128);

    // histograms storing median values based on the neighboring channels
    TH1F h_med_even("h_med_even", "h_med_even", 128, 0, 128);
    TH1F h_med_odd("h_med_odd", "h_med_odd", 128, 0, 128);

    // ---------------------------------------------- Histograms name
    h_ave_hits.SetTitle("Noise counts vs channel; Channel number; Counts");
    h_broken_ch.SetTitle("Noise counts vs channel; Channel number; Counts");
    h1_ave_even.SetTitle("Noise counts; Number of noise hits; Entries");
    h1_ave_odd.SetTitle("Noise counts; Number of noise hits; Entries");
    h2hits_amp.SetTitle("Noise counts distribution; Channel number; ADC value");
    h_thr_even.SetTitle("Threshold vs channel; Channel number; Threshold");
    h_thr_odd.SetTitle("Threshold vs channel; Channel number; Threshold");
    h_med_even.SetTitle(
        "Median of neighbors vs channel; Channel number; Median of neighbors");
    h_med_odd.SetTitle(
        "Median of neighbors vs channel; Channel number; Median of neighbors");

    int ch_hits[128] = {0};
    float ave_even = 0.;
    float ave_odd = 0.;

    // Initial iteration to calculate averages. Odd and Even channels separated
    for (int ch = 0; ch < 128; ch++) {
      for (int d = 1; d < 32; d++) {
        h2hits_amp.Fill(ch, 31 - d, vcnt[ch][d]);
        ch_hits[ch] += vcnt[ch][d];
      }
      h_ave_hits.Fill(ch, ch_hits[ch]);
      if (ch % 2 == 0) {
        h1_ave_even.Fill(ch_hits[ch]);
        h_ave_even.Fill(ch, ch_hits[ch]);
      } else {
        h1_ave_odd.Fill(ch_hits[ch]);
        h_ave_odd.Fill(ch, ch_hits[ch]);
      }
    }
    // medians of neigh_even, neigh_odd vectors
    float neigh_med_even, neigh_med_odd;
    // dynamic thresholds for even and odd channels based on the neighborhood
    float thr_even, thr_odd;

    // Iteration to discriminate broken channels
    for (int ch = 0; ch < 128; ch++) {

      // erase the neighbors of the previous channel
      neigh_odd.clear();
      neigh_even.clear();

      // the channel may not have enough neighbors to the left or right
      int start_ch = std::max(0, ch - n_neigh_ch);
      int end_ch = std::min(128, ch + n_neigh_ch + 1);

      if (ch % 2 == 0) {
        // if the channel is even

        for (int i = start_ch; i < end_ch; i++) {
          // Loop through neighboring channels to the left and right of the
          // current channel
          if (i != ch) {

            neigh_even.push_back(ch_hits[i]);
          }
        }

        // calculate the median of neighboring channels
        neigh_med_even = calculate_median(neigh_even, neigh_even.size());
        h_med_even.Fill(ch, neigh_med_even);
        thr_even = z_alpha * neigh_med_even;
        h_thr_even.Fill(ch, thr_even);

        if (ch_hits[ch] < thr_even) {
          // if the channel is broken based on the dynamic threshold

          broken_channels.push_back(ch);
          //______________________
          if (h_ave_hits.GetBinContent(ch + 1) == 0) {
            // No analog response (NAR): broken channel
            std::cout << "broken dynamic threshold NAR channel " << ch << '\n';
          } else {
            std::cout << "broken dynamic threshold channel " << ch << '\n';
          }
        }
      }

      else if (ch % 2 != 0) {
        // if the channel is odd

        for (int i = start_ch; i < end_ch; i++) {
          // Loop through neighboring channels to the left and right of the
          // current channel
          if (i != ch) {

            neigh_odd.push_back(ch_hits[i]);
          }
        }

        // calculate the median of neighboring channels
        neigh_med_odd = calculate_median(neigh_odd, neigh_odd.size());
        h_med_odd.Fill(ch, neigh_med_odd);
        thr_odd = z_alpha * neigh_med_odd;
        h_thr_odd.Fill(ch, thr_odd);

        if (ch_hits[ch] < thr_odd) {
          // if the channel is broken based on the dynamic threshold

          broken_channels.push_back(ch);
          if (h_ave_hits.GetBinContent(ch + 1) == 0) {
            // No analog response (NAR): broken channel
            std::cout << "broken dynamic threshold NAR channel " << ch << '\n';
          } else {
            std::cout << "broken dynamic threshold channel " << ch << '\n';
          }
        }
      }
    }

    // Writing the histogram of the broken, suspicious and noisy channels
    fill_hist_with_ch(&h_broken_ch, broken_channels, ch_hits);

    TCanvas *c1 = new TCanvas("c1", "c1", 900, 600);

    c1->SetGrid();
    h_ave_hits.Draw("HIST");
    h_ave_hits.SetStats(0);
    h_ave_hits.GetYaxis()->SetTitleOffset(1.2);
    h_ave_hits.GetYaxis()->SetNdivisions(505, 4);
    h_ave_hits.SetLineColor(kAzure - 3);
    h_ave_hits.GetYaxis()->SetRangeUser(0, 1.5 * h_ave_hits.GetMaximum());

    // Adding broken channels
    h_broken_ch.Draw("PHISTSAME");
    h_broken_ch.SetMarkerSize(1.5);
    h_broken_ch.SetMarkerStyle(23);
    h_broken_ch.SetMarkerColor(kBlack);

    h_med_even.Draw("PHISTSAME");
    h_med_even.SetLineWidth(0);
    h_med_even.SetMarkerStyle(kFullCircle);
    h_med_even.SetMarkerColor(kGreen + 2);

    h_thr_even.Draw("PHIST SAME");
    h_thr_even.SetLineWidth(0);
    h_thr_even.SetMarkerStyle(kFullDotLarge);
    h_thr_even.SetMarkerColor(kGreen + 2);

    h_med_odd.Draw("PHIST SAME");
    h_med_odd.SetLineWidth(0);
    h_med_odd.SetMarkerStyle(kFullCircle);
    h_med_odd.SetMarkerColor(kMagenta);

    h_thr_odd.Draw("PHIST   SAME");
    h_thr_odd.SetLineWidth(0);
    h_thr_odd.SetMarkerStyle(kFullDotLarge);
    h_thr_odd.SetMarkerColor(kMagenta + 2);

    TLegend *lg_brk_summ = new TLegend(0.65, 0.65, 0.90, 0.895, "", "nbNDC");
    lg_brk_summ->SetHeader("Connectivity check");
    lg_brk_summ->AddEntry(&h_ave_hits, "Noise hits", "l");
    lg_brk_summ->AddEntry(&h_broken_ch, "Broken channels", "p");
    lg_brk_summ->AddEntry(&h_med_even, "Median of neighboring even channels",
                          "p");
    lg_brk_summ->AddEntry(&h_med_odd, "Median of neighboring odd channels",
                          "p");
    lg_brk_summ->AddEntry(&h_thr_even, "Threshold broken even channels", "p");
    lg_brk_summ->AddEntry(&h_thr_odd, "Threshold broken odd channels", "p");

    lg_brk_summ->SetBorderSize(0);
    lg_brk_summ->SetFillColor(0);
    lg_brk_summ->Draw("");
    c1->Write();
    c1->Close();

    //  Writing histograms to file
    gStyle->SetOptStat();
    h_ave_hits.Write();
    h2hits_amp.Write();

    h_broken_ch.Write();

    h1_ave_even.Write();
    h1_ave_odd.Write();

    h_ave_even.Write();
    h_ave_odd.Write();

    // close the opened root file
    file1->Close();
    // Writing summary file
    write_results(filename_root, broken_channels);
    n_broken_channels += broken_channels.size();
    delete file1;
    delete c1;
    delete lg_brk_summ;
  }
  write_parameters(z_alpha, n_neigh_ch, n_broken_channels);
  return 0;
}
