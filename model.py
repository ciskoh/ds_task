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
    mic_to_soc: float = 0.5

    def calc(self, mic, doc, enz, soc, temp, **kwargs)->dict:
        """calculates DOC using Decomposition formula. WANRING: death is not calculated ans set equal to 0"""
        self.calc_v_max(temp)
        self.calc_assim(mic, doc)
        self.calc_decomp(enz, soc)
        return dict(self)
    
    def calc_v_max(self, temp):
        self.v_max = self.v_max_0 * exp(-self.ea /self.gas_const *(temp + 273))
    
    def calc_km(self, temp):
        self.km = self.km_slope * temp + self.km_0

    def calc_assim(self, mic, doc):
        self.assim = self.v_max * mic * (doc / (self.km + doc))
    
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
            - update vmax
            - update km_uptake
            - calc assim
        """
    
    def calc(self, temp, mic, doc, **kwargs):
        self.update_vmax_uptake(temp)
        self.update_km_uptake(temp)
        self.calc_assim(mic, doc)
        return dict(self)

    def update_vmax_uptake(self, temp):
        """calculates Arrhenius formula base don temperature and constants"""
        self.v_max_uptake = self.v_max_uptake_0 ^(- self.ea_uptake/self.gas_const * (temp + 273))
    
    def update_km_uptake(self, temp):
        self.km_uptake = self.km_uptake_slope * temp + self.km_uptake_0

    def update_assim(self, mic, doc ):
        self.assim = self.v_max_uptake * mic * (doc / (self.km_uptake + doc))

class CarbonUse(AbstractProcess):
    cue_slope: float = -0.016
    
    def calc(self, temp, cue, **kwargs):
        self.update_cue(temp, cue)
        return  dict(self)
    
    def update_cue(self, temp, cue):
        self.cue= self.cue_slope * temp + cue


class EnzymeProductionDecay(AbstractProcess):
    r_enz_prod: float = 0.000005
    r_enz_loss: float = 0.001
    
    def calc(self, mic, enz, **kwargs):
        self.calc_eprod(mic)
        self.calc_eloss(enz)
        return dict(self)
        
    def calc_eprod(self, mic):
        self.eprod = mic * self.r_enz_prod

    def calc_eloss(self, enz):
        self.eloss = enz * self.r_enz_loss


class MicDeath(AbstractProcess):
    r_death: float = 0.0002

    def calc(self, mic, **kwargs):
        self.death = self.r_death * mic
        return dict(self)


@dataclass
class ModelRunner:
    # processes
    encat: EnzymeCatalysis = field(default_factory=EnzymeCatalysis)
    uptake: Uptake = field(default_factory=Uptake)
    cbuse: CarbonUse = field(default_factory=CarbonUse)
    enzpd: EnzymeProductionDecay = field(default_factory=EnzymeProductionDecay)
    micdeath: MicDeath = field(default_factory=MicDeath)
    # system variables
    sysvar : dict = field(default_factory=dict)

    # extra intermediate variables, used for computation
    iv: dict =field(default_factory=dict)
    history: List[Dict[str, float]]= field(default_factory=list)

    def compute_model(self, dtemp, input_soc, input_doc):
        print(f"computing model for {dtemp=}")
        self.sysvar["temp"] += dtemp

        self.intermediate_calc(**self.sysvar)
        self.update_sysvar(input_soc, input_doc)

        print("computation finished")
        self.save_history()
        print(f"current state of the system:{self.sysvar=}")

    def intermediate_calc(self, **kwargs):
        self.iv={}
        for process in [self.encat,self.uptake, self.cbuse, self.enzpd, self.micdeath]:
            self.iv.update(process.calc(**kwargs)) 

    def update_sysvar(self, input_soc, input_doc):
        
        self.sysvar["soc"] += input_soc + self.iv["death"] * self.mic_to_soc - self.iv
        self.sysvar["doc"] += input_doc + self.iv["death"] * (1 - self.mic_to_soc) + self.iv["decomp"] + self.iv["eloss"] - self.iv["assim"]
        self.sysvar["mic"] += self.iv["assim"] * self.iv["cue"] - self.iv["death"] - self.iv["assim"]
        self.sysvar["enz"] += self.iv["eprod"] - self.iv["eloss"]
        self.sysvar["co2"] += self.iv["assim"] * (1-self.iv["cue"])

    def save_history(self):
        self.history.append(self.sysvar)

    def export_history(self, dst_path: str | p.Path):
        with open(dst_path, "w+") as f:
            json.dump(self.history, f)
            if p.Path(dst_path).is_file() and p.Path(dst_path).exists():
                print(f"computation history saved at {dst_path}")
            else:
                raise FileNotFoundError(
                    f"could not find history file at {dst_path}. Something went wrong"
                )
