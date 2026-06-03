// analyze_sec266.C — print histogram entries for secondary/absorber-2 diagnostics.
// Usage (from build/): root -l -b -q '../analyze_sec266.C("Hadr05.root")'
void analyze_sec266(const char* file = "Hadr05.root")
{
  TFile f(file);
  if (f.IsZombie()) {
    std::cerr << "Cannot open " << file << std::endl;
    return;
  }

  const char* partNames[] = {"alpha", "Li7", "Li6", "proton", "gamma", "e-", "C12"};

  auto printH1 = [&](const char* id, const char* label) {
    auto* h = dynamic_cast<TH1*>(f.Get(id));
    if (!h) {
      std::cout << "hist " << id << " (" << label << "): MISSING\n";
      return;
    }
    std::cout << "hist " << id << " (" << label << "): entries=" << h->GetEntries()
              << "  integral=" << h->Integral() << "  title=\"" << h->GetTitle() << "\"\n";
  };

  std::cout << "=== neutron interaction depth (absorber 2) ===\n";
  printH1("47", "elastic abs2");
  printH1("48", "inelastic+capture abs2");

  std::cout << "\n=== secondaries born in absorber 2 ===\n";
  printH1("266", "sec depth abs2");
  printH1("226", "capture gamma abs2");

  std::cout << "\n=== Edep by particle in absorber 2 (per-event fills; EndOfRun scaled) ===\n";
  for (int p = 0; p < 7; ++p) {
    printH1(Form("%d", 92 + p), partNames[p]);
  }
}
