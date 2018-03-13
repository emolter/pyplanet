def modify(gas, cloud, C, Cl):

    comment = "'Standard' Neptune atmospheric tweaking from Imke's code..."

    nAtm = len(gas[C['P']])
    for i in range(nAtm):

        # ## Process CO
        if gas[C['P']][i] > 0.1585:
            gas[C['CO']][i] = 0.0
        else:
            gas[C['CO']][i] = 1.0E-6
        # ## Process CO13
        gas[C['CO13']][i] = 1.0E-2 * gas[C['CO']][i]

        # ## Process HCN, SOLN, PH3
        gas[C['HCN']][i] = 0.0
        gas[C['SOLN']][i] = 0.0
        gas[C['PH3']][i] = 0.0

    return comment, gas, cloud
