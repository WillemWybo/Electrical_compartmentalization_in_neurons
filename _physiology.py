def setPhysiologyPas(phys_tree):
    # Major (2008) and Rhodes (2006) parameters
    for node in phys_tree:
        node.setPhysiology(0.8,      # Cm [uF/cm^2]
                           100./1e6, # Ra [MOhm*cm]
                          )
        node.addCurrent('L',  # leak current
                        100., # g_max [uS/cm^2]
                        -75., # e_rev [mV]
                       )

def setPhysiologyPasSpines(phys_tree):
    # Major (2008) and Rhodes (2006) parameters
    for node in phys_tree:
        if node.R < 0.6:
            # spine correction
            node.setPhysiology(0.8 * 1.92, # Cm [uF/cm^2]
                               100./1e6,   # Ra [MOhm*cm]
                              )
            node.addCurrent('L',         # leak current
                            100. * 1.92, # g_max [uS/cm^2]
                            -75.,        # e_rev [mV]
                           )
        else:
            node.setPhysiology(0.8,      # Cm [uF/cm^2]
                               100./1e6, # Ra [MOhm*cm]
                              )
            node.addCurrent('L',  # leak current
                            100., # g_max [uS/cm^2]
                            -75., # e_rev [mV]
                           )

def setPhysiologyActiveSoma(phys_tree):
    # Major (2008) and Rhodes (2006) parameters
    for node in phys_tree:
        node.setPhysiology(0.8,      # Cm [uF/cm^2]
                           100./1e6, # Ra [MOhm*cm]
                          )
        if node.index == 1:
            phys_tree[1].addCurrent('Kv3_1',  0.766     *1e6,   e_rev=-85., channel_storage=phys_tree.channel_storage)
            phys_tree[1].addCurrent('Na_Ta',  1.71      *1e6,   e_rev=50.,  channel_storage=phys_tree.channel_storage)
    phys_tree.fitLeakCurrent(e_eq_target=-75., tau_m_target=10.)