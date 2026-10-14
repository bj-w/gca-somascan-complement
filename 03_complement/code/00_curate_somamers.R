library(tidyverse)

res.limma <- read_delim("01_differential_abundance/data/limma.ref_dx.csv")

# lectin
lectin <- c(
  "MBL2_seq.3000.66",
  "FCN1_seq.22961.7",
  "FCN1_seq.3613.62",
  "FCN1_seq.15582.25",
  "FCN2_seq.3313.21",
  "FCN2_seq.13717.15",
  "FCN3_seq.14077.6",
  "FCN3_seq.5462.62",
  "MASP1_seq.8091.16",
  "MASP1_seq.3605.77",
  "MASP1_seq.9125.23"
)

# classical
classical <- c(
  "C1QC|C1QB|C1QA_seq.2753.2",
  "C1QC_seq.14100.63",
  "C1R_seq.3285.23",
  "C1S_seq.8840.61",
  "C1S_seq.20203.45"
)

# alternative
alternative <- c(
  "CFB_seq.4129.72",
  "CFD_seq.2946.52",
  "CFD_seq.13678.169",
  "CFP_seq.2960.66"
)

# C4
c4 <- c("C4B|C4A_seq.4481.34", "C4B|C4A_seq.18821.9", "C4B|C4A_seq.2182.54")

# C2
c2 <- c("C2_seq.3186.2")

# C3
c3 <- c(
  "C3_seq.2754.50",
  "C3_seq.4900.8",
  "C3_seq.2755.8",
  "C3_seq.4480.59",
  "C3_seq.2683.1",
  "C3_seq.5803.24"
)

# C5
c5 <- c("C5_seq.2381.52", "C5_seq.2851.63")

# membrane-attack complex
mac <- c(
  "C6|C5_seq.4482.66",
  "C6_seq.4127.75",
  "C7_seq.13731.14",
  "C7_seq.2888.49",
  "C8G|C8B|C8A_seq.2429.27",
  "C9_seq.3060.43",
  "C9_seq.13722.105"
)


# regulatory --------------------------------------------------------------
regulatory <- c(
  "SERPING1_seq.13710.6",
  "SERPING1_seq.4479.14",
  "C4BPA_seq.9449.150",
  "C4BPB_seq.7155.46",
  "CPN1_seq.7142.5",
  "CPN2_seq.6415.90",
  "CFHR5_seq.3666.17",
  "CFHR5_seq.16055.3",
  "CFHR5_seq.7885.17",
  "CFHR4_seq.6471.53",
  "CFHR1_seq.5982.50",
  "CFH_seq.4159.130",
  "CFHR3_seq.16836.1",
  "CFHR2_seq.15584.9",
  "CFHR1_seq.15468.14",
  "CFI_seq.2567.5",
  "VTN_seq.8280.238",
  "VTN_seq.13125.45",
  "VSIG4_seq.17678.28",
  "CR1_seq.19556.12",
  "CD55_seq.5069.9",
  "CD46_seq.17682.1",
  "CD59_seq.10480.33",
  "CD59_seq.11514.196",
  "CR2_seq.15570.99"
)

# combine and export ------------------------------------------------------
complement_soma <- bind_rows(
  data.frame(apt_id = lectin, section = "Lectin"),
  data.frame(apt_id = classical, section = "Classical"),
  data.frame(apt_id = alternative, section = "Alternative"),
  data.frame(apt_id = c4, section = "C4"),
  data.frame(apt_id = c2, section = "C2"),
  data.frame(apt_id = c3, section = "C3"),
  data.frame(apt_id = c5, section = "C5"),
  data.frame(apt_id = mac, section = "MAC"),
  data.frame(apt_id = regulatory, section = "Regulatory")
)
write_csv(complement_soma, "03_complement/data/complement_somamer_ids.csv")
