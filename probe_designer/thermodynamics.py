from Bio.SeqUtils import MeltingTemp as MT
import numpy as np


class ThermoParams:
    def __init__(self, temp=50.0, Na=50.0, Mg=0.0, dnac1=250.0, dnac2=250.0):
        """Class to set thermodynamic parameters of hybridization

        Parameters
        ----------
        temp : float, optional
            Temperature of reaction in Celsius
        Na : float, optional
            Monovalent cation concentration [mM]
        Mg : float, optional
            Divalent cation concentration [mM]
        dnac1 : float, optional
            Concentration of higher concentrated strand [nM]
        dnac2 : float, optional
            Concentration of lower concentrated strand [nM]
        """
        self.temp = temp
        self.Na = Na
        self.Mg = Mg
        self.dnac1 = dnac1
        self.dnac2 = dnac2
        self.paraG = self._create_deltaG_array(self.Na, self.temp)

    def __repr__(self):
        args = ', '.join(f'{k}={v}' for k, v in self.__dict__.items()
                         if k != 'paraG')
        return f'{self.__class__.__name__}({args})'

    def set(self, **kwargs):
        """Set thermodynamic attribute

        Parameters
        ----------
        **kwargs
            temp  (float): Temperature of reaction in Celsius
            Na    (float): Monovalent cation concentration [mM]
            Mg    (float): Divalent cation concentration [mM]
            dnac1 (float): Concentration of higher concentrated strand [nM]
            dnac2 (float): Concentration of lower concentrated strand [nM]

        Returns
        -------
        ThermoParams object
            Object containing thermodynamic parameters
        """
        keys_modified = kwargs.keys()
        for k, v in kwargs.items():
            try:
                v = float(v)
            except ValueError:
                raise ValueError(f'{k} expected float, received {type(v)}')
            setattr(self, k, v)
        if 'temp' in keys_modified or 'Na' in keys_modified:
            self.paraG = self._create_deltaG_array(self.Na, self.temp)
        return self

    def _create_deltaG_array(self, salinity, temp):
        """ Create paraG array for delta G calculations

        Parameters
        ----------
        salinity : float
            Monovalent cation concentration [mM]
        temp : float
            Temperature of reaction in Celsius

        Returns
        -------
        np.array
            Numpy array of paraG values [A T C G] x [A T C G]
        """
        salinity /= 1000
        temp += 273.15
        paraH = np.array([
            [-7.6, -7.2, -8.4, -7.8],
            [-7.2, -7.6, -8.2, -8.5],
            [-8.5, -7.8, -8, -10.6],
            [-8.2, -8.4, -9.8, -8]
        ])
        paraS = np.array([
            [-21.3, -20.4, -22.4, -21],
            [-21.3, -21.3, -22.2, -22.7],
            [-22.7, -21, -19.9, -27.2],
            [-22.2, -22.4, -24.4, -19.9],
        ])
        paraS = paraS + 0.368 * np.log(salinity)
        paraG = paraH - temp * paraS / 1000
        return paraG


thermo = ThermoParams()


def calc_deltaG(sequence, temp=thermo.temp, paraG=thermo.paraG):
    """Calculate deltaG from sequence. Default paraG is calculated
    using the ThermoParams default values. If a different temperature/salinity
    is required, create a new ThermoParams object and use class temp and paraG
    as arguments

    Parameters
    ----------
    sequence : str
        Sequence of ACTG bases
    temp : float, optional
        Temperature of reaction in Celsius
    paraG : list of list, optional
        ParaG array

    Returns
    -------
    float
        Delta G estimation of sequence (or np.nan if fails due to
        non-ACTG base in sequence)
    """
    temp += 273.15
    # Initialization energy
    dG = 0.2 - temp * -5.7 / 1000

    mapping = {'A': 0, 'T': 1, 'C': 2, 'G': 3}
    for base1, base2 in zip(sequence, sequence[1:]):
        a = mapping.get(base1.upper(), None)
        b = mapping.get(base2.upper(), None)
        if a is None or b is None:
            return np.nan
        dG += paraG[a, b]
    return dG


def calc_tm(sequence, dnac1=250, dnac2=250, Na=50, Mg=0, dNTPs=0, c_seq=None):
    """ Return melting temperture as float or np.nan if fails
    (thermo not available)

    Parameters
    ----------
    Na : float, optional
        Monovalent cation concentration [mM]
    Mg : float, optional
        Divalent cation concentration [mM]
    dnac1 : float, optional
        Concentration of higher concentrated strand [nM]
    dnac2 : float, optional
        Concentration of lower concentrated strand [nM]
    dNTPs : float, optional
        Concentration of dNTPs [mM]
     """
    try:
        return MT.Tm_NN(sequence, dnac1=dnac1, dnac2=dnac2, Na=Na,
                        Mg=Mg, dNTPs=dNTPs, c_seq=c_seq)
    except ValueError:
        return np.nan


def calc_gc(sequence):
    """ Return float of GC content """
    return sum(x in 'GC' for x in sequence) / len(sequence)
