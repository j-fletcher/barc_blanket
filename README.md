# barc_blanket


## Install environment

```
conda env create -f environment.yml
```

```
conda activate barc-blanket-env
```

## Configure cross sections

Certain tank analytes are not found in standard neutron cross section libraries. As a temporary measure, we will tell OpenMC to use TENDL-2019 cross sections for these analytes, alongside the default library (in our case, ENDF/B-VIII.0). It is advisable that users make a note for future projects that we have modified our installation of the standard cross section data! First, locate your cross_sections.xml:

```
echo $OPENMC_CROSS_SECTIONS
```

Finding this file will also show you where OpenMC is being directed to look for cross sections. There should be a "neutron" folder in its parent directory (if entries in cross_sections.xml are formatted as described below); this is where we will put the TENDL cross sections.

The new cross sections can be found on the repo under TENDL_cross_sections. Copy these into the "neutron" folder from above:

```
Ba137_m1.h5
C14.h5
Cd113_m1.h5
Co60.h5
Nb93_m1.h5
Ra228.h5
```

Next, within the <cross_sections> section of cross_sections.xml, add the following lines:

```
<library materials="Cd113_m1" path="neutron/Cd113_m1.h5" type="neutron" />
<library materials="Ba137_m1" path="neutron/Ba137_m1.h5" type="neutron" />
<library materials="C14" path="neutron/C14.h5" type="neutron" />
<library materials="Co60" path="neutron/Co60.h5" type="neutron" />
<library materials="Nb93_m1" path="neutron/Nb93_m1.h5" type="neutron" />
<library materials="Ra228" path="neutron/Ra228.h5" type="neutron" />
```

Save your progress. OpenMC should now be able to handle the more "exotic" waste products.
