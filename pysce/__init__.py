####################################################################################################
# # Copyright (C) 2024-Present - Daniel Charytonowicz - All Rights Reserved
# # Written by: Daniel Charytonowicz
# # Contact: daniel.charytonowicz@icahn.mssm.edu
# ###################################################################################################

from ._velocity import score_angular_velocity_entropy, ensure_velocity


def __getattr__(name):
    # Expression entropy pulls in torch/scanpy/networkx; import it lazily so
    # the velocity tools work without those installed.
    if name == 'score':
        from ._tools import score
        return score
    raise AttributeError(f"module 'pysce' has no attribute {name!r}")


__all__ = ['score', 'score_angular_velocity_entropy', 'ensure_velocity']
