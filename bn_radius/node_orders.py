"""Original (BIF-file) variable order of the bnlearn networks.

The benchmark query sampler of the paper draws vertex *indices*, so the results depend on this order.
Storing it explicitly makes the sampled results independent of the pgmpy version."""

BIF_ORDER = {
    "asia": [
        'asia', 'tub', 'smoke', 'lung', 'bronc', 'either',
        'xray', 'dysp',
    ],
    "insurance": [
        'GoodStudent', 'Age', 'SocioEcon', 'RiskAversion', 'VehicleYear', 'ThisCarDam',
        'RuggedAuto', 'Accident', 'MakeModel', 'DrivQuality', 'Mileage', 'Antilock',
        'DrivingSkill', 'SeniorTrain', 'ThisCarCost', 'Theft', 'CarValue', 'HomeBase',
        'AntiTheft', 'PropCost', 'OtherCarCost', 'OtherCar', 'MedCost', 'Cushioning',
        'Airbag', 'ILiCost', 'DrivHist',
    ],
    "alarm": [
        'HISTORY', 'CVP', 'PCWP', 'HYPOVOLEMIA', 'LVEDVOLUME', 'LVFAILURE',
        'STROKEVOLUME', 'ERRLOWOUTPUT', 'HRBP', 'HREKG', 'ERRCAUTER', 'HRSAT',
        'INSUFFANESTH', 'ANAPHYLAXIS', 'TPR', 'EXPCO2', 'KINKEDTUBE', 'MINVOL',
        'FIO2', 'PVSAT', 'SAO2', 'PAP', 'PULMEMBOLUS', 'SHUNT',
        'INTUBATION', 'PRESS', 'DISCONNECT', 'MINVOLSET', 'VENTMACH', 'VENTTUBE',
        'VENTLUNG', 'VENTALV', 'ARTCO2', 'CATECHOL', 'HR', 'CO',
        'BP',
    ],
}
