""" 
    This contains the classes 
    that compose the Enzyme based model 
    presented by Allison et al. 201067
"""
from abc import ABC, abstractmethod
from dataclasses import dataclass, field
from typing import Dict, List
import pathlib as p
from collections import Counter
import json
from math import exp
from decimal import Decimal, getcontext


@dataclass
class AbstractProcess(ABC):

    @abstractmethod
    def calc(self, input):
        pass


class EnzymeCatalysis(AbstractProcess):
    v_max_0: int = 100000000
    ea= 47
    gas_const = 0.008314
    km_0: int = 500.0
    km_slope = 5

    def calc(self,  enz, soc, temp, **kwargs)->dict:
        """calculates DOC using Decomposition formula. WANRING: death is not calculated ans set equal to 0"""
        self.calc_v_max(temp)
        self.calc_km(temp)
        self.calc_decomp(enz, soc)
        return self.__dict__
    
    
    def calc_v_max(self, temp):
        getcontext().prec = 28 # set decimal precision

        fract = Decimal(str(self.ea)) / (Decimal(self.gas_const) * (Decimal(str(temp)) + Decimal('273')))
        exp_factor = Decimal(str(exp(-fract)))
        self.v_max = float(Decimal(str(self.v_max_0)) * exp_factor)
    
    def calc_km(self, temp):
        self.km = self.km_slope * temp + self.km_0
    
    def calc_decomp(self, enz, soc):
        self.decomp = self.v_max * enz * (soc / (self.km + soc))


class Uptake(AbstractProcess):
    v_max_uptake_0: int = 100000000
    km_uptake_0: float = 0.1
    ea_uptake: float = 47
    gas_const: float = 0.008314
    km_uptake_slope: float = 0.01


    def calc(self, temp, mic, doc, cue, death, eprod, **kwargs):
        """calculates dMIC using microbial mass change formula
        steps:
            - calc vmax_uptake
            - calc km_uptake
            - calc assim
        """
    
    def calc(self, temp, mic, doc, **kwargs):
        self.calc_vmax_uptake(temp)
        self.calc_km_uptake(temp)
        self.calc_assim(mic, doc)
        return self.__dict__

    def calc_vmax_uptake(self, temp):
        """calculates Arrhenius formula based 
        on temperature and constants
        Uses decimal library to account for very small exponential
        """

        fract = Decimal(str(self.ea_uptake)) / (Decimal(self.gas_const) * (Decimal(str(temp)) + Decimal('273')))
        exp_factor = Decimal(str(exp(-fract)))
        self.v_max_uptake = float(Decimal(str(self.v_max_uptake_0)) * exp_factor)
    

    def calc_km_uptake(self, temp):
        self.km_uptake = self.km_uptake_slope * temp + self.km_uptake_0

    def calc_assim(self, mic, doc ):
        self.assim = self.v_max_uptake * mic * (doc / (self.km_uptake + doc))

class CarbonUse(AbstractProcess):
    cue_slope: float = -0.016
    cue_0: float=0.63
    
    def calc(self, temp, **kwargs):
        self.calc_cue(temp)
        return self.__dict__
    
    def calc_cue(self, temp):
        self.cue= self.cue_slope * temp + self.cue_0


class EnzymeProductionDecay(AbstractProcess):
    r_enz_prod: float = 0.000005
    r_enz_loss: float = 0.001
    
    def calc(self, mic, enz, **kwargs):
        self.calc_eprod(mic)
        self.calc_eloss(enz)
        return self.__dict__
        
    def calc_eprod(self, mic):
        self.eprod = mic * self.r_enz_prod

    def calc_eloss(self, enz):
        self.eloss = enz * self.r_enz_loss


class MicDeath(AbstractProcess):
    r_death: float = 0.0002

    def calc(self, mic, **kwargs):
        self.death = self.r_death * mic
        return self.__dict__


@dataclass
class ModelRunner:
    
    # system variables define the inital /current state of the system are exported to history
    sysvar : dict = field(default_factory=dict)
    # parameters for SOC /DOC /ENZ / CO2  recalc
    mic_to_soc: float=0.5
    # processes
    encat: EnzymeCatalysis = field(default_factory=EnzymeCatalysis)
    uptake: Uptake = field(default_factory=Uptake)
    cbuse: CarbonUse = field(default_factory=CarbonUse)
    enzpd: EnzymeProductionDecay = field(default_factory=EnzymeProductionDecay)
    micdeath: MicDeath = field(default_factory=MicDeath)

    # extra intermediate variables, used for computation
    iv: dict =field(default_factory=dict)
    history: List[Dict[str, float]]= field(default_factory=list)

    def compute_model(self, dtemp, input_soc, input_doc, verbose=True):
        if verbose:
            print(f"computing model for {dtemp=}, {input_soc=},{input_doc=} ")
            print(f"inital system state: {self.sysvar=}")
        self.sysvar["temp"] += dtemp

        self.calc_subprocesses(**self.sysvar)
        self.update_sysvar(input_soc, input_doc)

        self.save_history()
        if verbose:
            print("computation finished")
            print(f"current state of the system:{self.sysvar=}")

    def calc_subprocesses(self, **kwargs):
        self.iv={}  
        for process in [self.encat,self.uptake, self.cbuse, self.enzpd, self.micdeath]:
            self.iv.update(process.calc(**kwargs)) 

    def update_sysvar(self, input_soc, input_doc):
        
        self.sysvar["soc"] += input_soc + self.iv["death"] * self.mic_to_soc - self.iv["decomp"]
        self.sysvar["doc"] += input_doc + self.iv["death"] * (1 - self.mic_to_soc) + self.iv["decomp"] + self.iv["eloss"] - self.iv["assim"]
        self.sysvar["mic"] += self.iv["assim"] * self.iv["cue"] - self.iv["death"] - self.iv["eprod"]
        self.sysvar["enz"] += self.iv["eprod"] - self.iv["eloss"]
        self.sysvar["co2"] += self.iv["assim"] * (1-self.iv["cue"])

    def save_history(self):
        self.history.append(self.sysvar.copy())

    def export_history(self, dst_path: str | p.Path):
        with open(dst_path, "w+") as f:
            json.dump(self.history, f)
            if p.Path(dst_path).is_file() and p.Path(dst_path).exists():
                print(f"computation history saved at {dst_path}")
            else:
                raise FileNotFoundError(
                    f"could not find history file at {dst_path}. Something went wrong"
                )
