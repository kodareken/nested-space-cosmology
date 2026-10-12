"""Independent post-result binder for the completed TDG10 QA1 terminal.

The one-time live path authenticates the immutable QA1 authority commit and
every committed FRZ1 delta blob, reads the canonical nested four-leaf raw
namespace with no-following race-checked opens, restores generation 9 /
journal sequence 10 / retry 3 from low-level primitives, constructs the
seven SSPRK3-on-inherited-SBP4 shadows, and recomputes all eighteen exact
complete-C channels with the frozen dual-route owner.  It deliberately does
not import the QA1 runner, the TI2 runner, or QA1 authority decisions, and
it never restores generation 10.  Ordinary verification consumes only the
tracked compact certificate and is raw/store/shadow blind.
"""

from __future__ import annotations

from copy import deepcopy
from dataclasses import asdict, dataclass, is_dataclass
from fractions import Fraction
from hashlib import sha256
from io import BytesIO
import json
import math
import os
from pathlib import Path
import stat
import subprocess
import tomllib
from typing import Any, Mapping, NoReturn, Sequence
import zipfile

import numpy as np

from . import hlt16_campaign_runtime as campaign_runtime
from . import proto15_runtime as p15
from . import tdg6_temporal_admission_runtime as tdg6
from .hlt16_campaign_store import HLT16CampaignStore
from .hlt16_member_codec import ARRAY_NAMES
from .numerical_engine import COMPARATOR_METHOD, PRIMARY_METHOD, array_content_sha256
from .proto19_gr0_static_factory import build_static_gr0_shells
from .tdg6_temporal_admission_design import (
    CertifiedMagnitudeInterval,
    TDG6_COMPLETE_STATE_CHANNELS,
    TDG6_ORDER_SQUARED_MULTIPLIER,
    classify_tdg6_channel,
    three_halves_order_passes_squared,
)
from .tdg9_exact_temporal_arithmetic import (
    exact_hermite_coefficients,
    restrict_exact_cubic_to_half,
    subtract_exact_cubics,
)
from .tdg9_local_extrema import (
    EVALUATOR_ID as PRIMARY_LOCALIZER_ID,
    LocalCubic,
    localize_absolute_maximum,
)
from .tdg9_local_extrema import _stationary_intervals as _primary_stationary_intervals
from .tdg9_local_extrema_independent_v2 import (
    EVALUATOR_ID as INDEPENDENT_LOCALIZER_ID,
    IndependentLocalCubicV2,
    RootIsolationInconclusive,
    localize_absolute_maximum_independently_v2,
)


SCHEMA_VERSION = 1
ARTIFACT_ID = "FGC-1-TDG10-QA1-PREF1"
CLASSIFICATION = (
    "independently_bound_retry3_ssprk3_sbp4_exact_complete_C_all_channel_pass_terminal"
)
CONFIG_PATH = "configs/fgc/fgc-1-tdg10-qa1-pref1.toml"
RESULT_PATH = "results/fgc-1-tdg10-qa1-pref1.json"
OWNER_DOCUMENT = "docs/fgc-tdg10-qa1-pref1.md"
TARGET_PROTOCOL = "FGC-2-SF1-PROTO18"
COMPACT_RESULT_SHA256 = (
    "3a107bcede479df4326aac5e194beefa1042860a8be498d7dc5d0519fdbc8603"
)

AUTHORITY_COMMIT = "aeaf0498fd37a51d0e2c7efcf695147f7a9822e7"
AUTHORITY_PARENT = "49c514ac283ae3ec8091a6638442da6fdaf58db0"
RAW_NAMESPACE = "runs/fgc-2-sf1/tdg10-qa1/retry3-ssprk3-sbp4-exact-complete-c"
RAW_SCHEMA = "FGC-1-TDG10-QA1-raw-v1"
RAW_CLASSIFICATION = "completed_exact_all_channel_pass_no_state_advance"
NONPASS_CLASS = "completed_exact_all_channel_nonpass_no_state_advance"
RUNNER_ID = "FGC-1-TDG10-QA1-RUN1"
QA1_ARTIFACT_ID = "FGC-1-TDG10-QA1-FRZ1"
DIAGNOSTIC_ENDPOINT_SCHEMA = "FGC-1-TDG10-QA1-diagnostic-fine-endpoint-v1"
DIAGNOSTIC_CLASSIFICATION = "diagnostic_nonaccepted_complete_fine_endpoint"
RAW_MANIFEST_SHA256 = "cdea5f1148486b6d0c2c5a3fee222a224ede758680e481efb8d6a78672c52b36"
RAW_TERMINAL_SHA256 = "c0aea714469dd79874220f12011fb6c46a8468c7d5866fcffa6daffe06427c98"
DESCRIPTOR_RAW_SHA256 = (
    "2c4a560c87d08f43dec5644d5af65c44663fbd77cb07476742f58a9a0e448532"
)
DESCRIPTOR_CONTENT_ADDRESS = (
    "69478ba14677cec5dfcc7a2f802f2e5a5fe729907a225a42b560d489c55e2c3e"
)
PAYLOAD_RAW_SHA256 = "740b371a9720b53e28164d67d22364172df17a2a30a738043c12033a722eed4f"
PAYLOAD_SEMANTIC_SHA256 = (
    "ea961ef96a03ccfd572c53e5bfac56cd9898cbe464a5f10d928022c9a7359a9e"
)
DESCRIPTOR_RELATIVE = f"fine-endpoint/{DESCRIPTOR_CONTENT_ADDRESS}.json"
PAYLOAD_RELATIVE = f"payloads/{PAYLOAD_SEMANTIC_SHA256}.npz"

STORE_PATH = "runs/fgc-2-sf1/tdg8-rcv3/calibration"
SEALED_STORE_LEAF_COUNT = 115
SEALED_STORE_SHA256 = "5917d70de3080dcc6b7abeb7fdb7cc3370c6140cbeb37bde0c142d6a2d36a445"
UR1_PREF1_RESULT_SHA256 = (
    "27c2aef9b033092285b868f583441842de99f61941ff57aa7e62749cc8d91300"
)
EXPERIMENT_LABEL = (
    "retry3_SSPRK3_on_inherited_SBP4_exact_complete_C_all_18_channel_qualification"
)
TABLEAU_RUNTIME_SELECTOR = COMPARATOR_METHOD
INTERVAL_OWNER = "exact_radius_free_complete_C_dual_rational_localizer"
MEMBER_KEY = "RK4-2049"
PHYSICAL_STATE_SHA256 = (
    "3cd9f576d079595343e8cdaabd33dd0153763b88d1beecd7d85717bdba4b6c9a"
)
ACCEPTED_TIME_HEX = "0x1.78554de5a30e0p+0"
POINT_COUNT = 2049
OWNED_ROW_COUNT = 2044
RETRY = 3
PRIOR_RETRY_COUNT = 2
WIDTH_HEX = "0x1.aaa9612df8000p-11"
MEMBER_DESCRIPTOR_SHA256 = (
    "77847126340e78c8ac300fac2795bceda724c15e0894e34444b0da42a84166b4"
)
COORDINATES_SHA256 = "2a347b314de5e1e94ac0c999db283e03498070ece6762c53252fa41c1e04a646"
GRID_SPACING_HEX = "0x1.0000000000000p-4"
OUTER_RADIUS_HEX = "0x1.0000000000000p+7"
TRANSACTION_SHA256 = "7a90163d2eb4252fa7a1bbdf55d12f92a9127ed7cba5a37c28456e72d732b22b"
HISTORICAL_JOURNAL_SHA256 = (
    "0d21f650ccc46db778fd81fbebf5acab6394e1c2f7972b7940e24c76af9e1f59"
)
PREDECESSOR_GENERATION = 9
FORBIDDEN_RETRY3_GENERATION = 10
CHECKPOINT_SHA256 = "eb6fddc480c94aa7ed15c80399fc2b26637693efef02399f267075fc6c258e56"
CHECKPOINT_RAW_SHA256 = (
    "07572d4450829ff2e2bad65e217a917b0d31ba677a3b8fc20a3be43e849ba846"
)
CHECKPOINT_JOURNAL_SEQUENCE = 10
CHECKPOINT_JOURNAL_SHA256 = (
    "5b533eb7009a9c7c3f353d26813f9cbbde8574afae77492929b435b4c701c5ff"
)
CHECKPOINT_JOURNAL_RAW_SHA256 = (
    "b50dc39ebd8bbdd3729d40d9e9ae4b223ff2ef19389d4ca74f9a316bac72202f"
)
HISTORICAL_RETRY3_REJECTION_SEQUENCE = 11
HISTORICAL_RETRY3_REJECTION_SHA256 = (
    "341cd8cd328434d85774bcb452889ac1886c520a72e945a92337b72cf564161a"
)
HISTORICAL_RETRY3_REJECTION_RAW_SHA256 = (
    "e4f9ca42f2d014d4bdbd875a0f24a9bbb9c04159fb0d9e0e7de8d0d72db61fdf"
)
PRIMARY_REFINEMENT_DEPTH = 160
PER_CHANNEL_MAXIMUM_CANDIDATES_D01 = 16352
PER_CHANNEL_MAXIMUM_CANDIDATES_D12 = 32704
D01_POLYNOMIAL_COUNT = 2 * OWNED_ROW_COUNT
D12_POLYNOMIAL_COUNT = 4 * OWNED_ROW_COUNT
CHANNEL_ORDER = TDG6_COMPLETE_STATE_CHANNELS
SUFFICIENT_CONDITION = "8*U12^2 <= L01^2"
PRIMARY_EVALUATOR_ID = "tdg9_loc1_derivative_monotone_bisection_v1"
INDEPENDENT_EVALUATOR_ID = "tdg9_loc2_endpoint_deflated_adaptive_discriminant_v2"
TDG10_EVALUATOR_ID = "tdg10_exact_complete_c_dual_route_v1"
_ROW_HASH_DOMAIN = b"TDG10-EXACT-COMPLETE-C-ROW-STREAM-v1\n"
_COEFFICIENT_HASH_DOMAIN = b"TDG10-EXACT-COMPLETE-C-CUBIC-STREAM-v1\n"
_COMBINED_HASH_DOMAIN = b"TDG10-EXACT-COMPLETE-C-COMBINED-v1\n"
_SURVIVOR_KEY_HASH_DOMAIN = b"TDG10-EXACT-COMPLETE-C-SURVIVOR-KEY-STREAM-v1\n"
_STATIONARY_COUNT_DOMAIN = b"TDG9-LOC2-STATIONARY-COUNT-STREAM-v1\n"
RETRY3_REPLAY = {
    "retry": RETRY,
    "prior_retry_count": PRIOR_RETRY_COUNT,
    "predecessor_generation": PREDECESSOR_GENERATION,
    "checkpoint_sha256": CHECKPOINT_SHA256,
    "checkpoint_raw_sha256": CHECKPOINT_RAW_SHA256,
    "journal_sequence": HISTORICAL_RETRY3_REJECTION_SEQUENCE,
    "journal_sha256": HISTORICAL_RETRY3_REJECTION_SHA256,
    "journal_raw_sha256": HISTORICAL_RETRY3_REJECTION_RAW_SHA256,
    "attempted_width_hex": WIDTH_HEX,
    "transaction_sha256": TRANSACTION_SHA256,
}
WORK_BUDGET = {
    "SSPRK3_records_per_proposal": 4,
    "channel_count": 18,
    "fourth_width_authorized": False,
    "generation_10_authorized": False,
    "maximum_stage_and_endpoint_RHS_records": 28,
    "owned_row_count": OWNED_ROW_COUNT,
    "per_channel_maximum_candidates_D01": PER_CHANNEL_MAXIMUM_CANDIDATES_D01,
    "per_channel_maximum_candidates_D12": PER_CHANNEL_MAXIMUM_CANDIDATES_D12,
    "primary_refinement_depth": PRIMARY_REFINEMENT_DEPTH,
    "resource_escalation_authorized": False,
    "retry": 3,
    "retry_4_authorized": False,
    "retry_5_authorized": False,
    "retry_count": 1,
    "shadow_paths": 7,
    "shadow_proposals": 7,
}
EXPECTED_REPLAY_RECEIPT = {
    "accepted_time_hex": ACCEPTED_TIME_HEX,
    "attempted_width_hex": WIDTH_HEX,
    "checkpoint_journal_sequence": CHECKPOINT_JOURNAL_SEQUENCE,
    "descriptor_sha256": MEMBER_DESCRIPTOR_SHA256,
    "historical_journal_sha256": HISTORICAL_JOURNAL_SHA256,
    "historical_retry3_rejection_sequence": (HISTORICAL_RETRY3_REJECTION_SEQUENCE),
    "member_key": MEMBER_KEY,
    "predecessor_generation": PREDECESSOR_GENERATION,
    "retry": RETRY,
    "state_sha256": PHYSICAL_STATE_SHA256,
    "transaction_sha256": TRANSACTION_SHA256,
}
EXPECTED_OUTPUT_LEAVES = (
    "manifest.json",
    "terminal.json",
    DESCRIPTOR_RELATIVE,
    PAYLOAD_RELATIVE,
)
EXPECTED_D12_UPPERS: tuple[tuple[str, str, str], ...] = (
    ("u:alpha", "7", "4503599627370496"),
    (
        "u:v",
        "831246704379425320646844156357345306227378210667277802810441586525544667119564294384172239380933057464118186088118406546553196682415358116335925762331592405796283716679917242735941978841078851146072504810714778883829204990073",
        "295750865823626734935838402260649096162709258871437468708140147072020998305800988936042902977479998441890691021854653331929270587495459189600629109461312184486386857586641474561234081783376970105371385333130499991410214954938786948239589376",
    ),
    (
        "u:lambda",
        "941904875588187868609460366327497949502423693677001520645788018745683247265584348955205512863477656644866626540443910769233099251018320754553856524344196212143961548380533780173922514799056094335310049327113077499406402813761",
        "302996581447774369040660351219870748373769124340360277768358163275803692410361074457916049659363433165411007232459726777902855325939042939604327990711064915108572756450512843412031207265601739588570594294928279595541073873072181037762609152",
    ),
    (
        "u:R",
        "90573206688753319135650781751004232376815226772476669667095",
        "1556021749743481858702452193862444089152448363676305536488467735152427008",
    ),
    ("u:phi", "42561", "9444732965739290427392"),
    ("u:chi", "36923", "144115188075855872"),
    (
        "p:alpha",
        "738509182912386578274345913129263367120556497067197951058926583656939543073370986246419167950418325317249402315952329019663891171113225924157388706317992755722407117843303209321401638448253751472815983157634460394498812372902658037",
        "43107138964873214920915984488608258936653938829903724629129002783481806955641076626111465123895827284872313181779077956090379599272156658449229231148704167706035027637381265666507818945011490296125356699723163176230836110831740512003378504531968",
    ),
    ("p:v", "30501", "1152921504606846976"),
    (
        "p:lambda",
        "7014140437795328040402718299826859973641382862792043150123048748647173582096232799599011707256591426844833571483304291770111544198744274548333807369142136694009687929371150751432178057472501078247077261673044704609450405",
        "267666884339797630337074134275806994718212944386620236616931701614301192126851464704197973022679782872446119884526704366699825083594329790230098960885049437749894363088017880700848599332271103063752428018765660469274847514425196281856",
    ),
    (
        "p:R",
        "205453089205130792772960522267762581495423099731902652844690509738572450588064272781441534428348340774325183247038626695113730125049186883175679724328559168541161656508008220181279103780208205676084367248445422045019807764194022041",
        "1020878746438771740852398381781227899324268273400524045226350631716926845118843567519452196065909963640519615911191170690162467952877250205474468436448970757581442014131881697004844521751500807682205614898680851119108040239967528437602347646976",
    ),
    (
        "p:phi",
        "32020759961110031496184202084270718702590309696484250145492625721265131609234182854958127224187421638387678320054460109021096445329322964767126856200918181704925630003650744420703224516159075313453979378813171132337813001571713",
        "635662448629622990806484057378795784400552398515919715004861941090940239166619966218965800235756536906530535740668941298420636436851196018271445719060580799752574569857402470798106188499568666041069767046314444804542019792030091279851172397056",
    ),
    ("p:chi", "35769", "9007199254740992"),
    (
        "q:alpha",
        "677671153093073804992540777508262512056066988643646584919699038822474085214166021934679846080113276624179501266459049393346882764141956828596852324518270046688023469258579427184421076944771855545566681296970775276254504375400121247",
        "49686449661782400195837265938457615647147731320628329162112268129051730074186585942732679450387594150867988460410910418762726745486612512176284323311978823072763140661683893589700255389308725589460273144640869219698582792013904061314654898487296",
    ),
    (
        "q:v",
        "100149261937569985961112547373601025669002092579037065162063314884002370320178028943967396319501574717506893186394278769140962549639397710172983919270327228644991667587882288396868816990124416201825251662845199023407814295063204563",
        "2791330633403527773992201078929220521872144201802688962882261076179167612484515146938144284250987425457925315106188391955304001433864172297419534260288153737384486676555273148762832280809918960828755445532413167673422379474624594394189555302400",
    ),
    (
        "q:lambda",
        "2281590850079873425396997939245464842298933988105131072176133744991992473873642140950656049858696009496454535229813677567919213842442515068833828905243488171488226505094273606086220487073644019952690319551567015436408027441190523",
        "36196799313590694054591706950528377265915990160225066261019437943587203358127317806198935948710079739786691272986622378759196979181300820545971266116788573957273722996954494620401610929354988628027717033933675548188978497047288620379013120000",
    ),
    ("q:R", "813", "4503599627370496"),
    (
        "q:phi",
        "399872800823967379376145613340427973136718907265319098088552504400065303991287306648318245272754005584903408351883610086670119024893800990086937150737004001227636538279838697532305515154563374397471871693202551147325903310302311",
        "8131951761954608862097558052099019851800807105663855319559682883392099109736665793427856502444966677781914050905880936166282861776671894035260963984027731118724315998576235170328830015664225156205175749381740674715201202934852325779268287070208",
    ),
    (
        "q:chi",
        "374095914125166599532290292242175031337597454840220267608729866811855725505131318443738011560152237760320471538067578258302765581112290218344659851612388770090011476306631362555524891489712049401975883406264218012823151236924657",
        "94147985503552033727312975695063760296660479504355391268057087117324167527805167620463724987347067596421860477943928463984315947944132750586098152369184049547396002130790033574480928471216364280774158302728170917953774931539217825048559616",
    ),
)
EXPECTED_PAYLOAD_MANIFEST = (
    {
        "byte_count": 98352,
        "bytes_sha256": (
            "d19cf412453c40a4fc652d2688bf8bc134428cf943a2d04805f9ca73f1f07d04"
        ),
        "dtype": "<f8",
        "name": "u",
        "order": "C",
        "shape": [2049, 6],
    },
    {
        "byte_count": 98352,
        "bytes_sha256": (
            "5390180ef14855aa4ceda1cbfaf5eb5a443c730d63f9aaccb78f8e0a636352aa"
        ),
        "dtype": "<f8",
        "name": "p",
        "order": "C",
        "shape": [2049, 6],
    },
    {
        "byte_count": 98352,
        "bytes_sha256": (
            "ab0bdd8f31bdcc64fd6674788be11f8646b85079b38e07341069f5daf79804a7"
        ),
        "dtype": "<f8",
        "name": "q",
        "order": "C",
        "shape": [2049, 6],
    },
    {
        "byte_count": 16392,
        "bytes_sha256": (
            "7bc5799c10eaefc1ba90acac2fd238fbefec192cc92ba269694f96ee48d26626"
        ),
        "dtype": "<f8",
        "name": "grid_coordinates",
        "order": "C",
        "shape": [2049],
    },
    {
        "byte_count": 384,
        "bytes_sha256": (
            "3b06bbb6620f366a5e53886a3d20b05d9c97d21f94e99c17d471daf1f94898ae"
        ),
        "dtype": "<f8",
        "name": "tracer_labels",
        "order": "C",
        "shape": [48],
    },
    {
        "byte_count": 384,
        "bytes_sha256": (
            "d81ab9b5f3925ad9f14604ac7dc69acab7f44d1c763332c1f3a2f6f292247a7a"
        ),
        "dtype": "<f8",
        "name": "tracer_positions",
        "order": "C",
        "shape": [48],
    },
    {
        "byte_count": 384,
        "bytes_sha256": (
            "46e5faab73ca92e8fd06ec2ac6ae965deacc93d532f728e5e19dfab3d6794993"
        ),
        "dtype": "<f8",
        "name": "tracer_proper_times",
        "order": "C",
        "shape": [48],
    },
    {
        "byte_count": 9216,
        "bytes_sha256": (
            "a614094990fc56c8436dd47f31f6df20e050eb23cd9ae1d009caed3fd9cfdbef"
        ),
        "dtype": "<f8",
        "name": "event_proper_times",
        "order": "C",
        "shape": [24, 48],
    },
    {
        "byte_count": 55296,
        "bytes_sha256": (
            "cc0161997901dc786289370c4ce1c33945eb781eb04a367f07a53d3113140064"
        ),
        "dtype": "<f8",
        "name": "event_fields",
        "order": "C",
        "shape": [24, 48, 6],
    },
)
ENDPOINT_FALSE_FLAGS = (
    "FGCQR_authorized",
    "GR0_calibration_authorized",
    "GR0_calibration_completed",
    "PDE_state_commit_authorized",
    "PDE_state_committed",
    "SGBL_authorized",
    "accepted_state",
    "campaign_descriptor",
    "campaign_state_write_authorized",
    "candidate_branch_opened",
    "candidate_execution_authorized",
    "checkpoint",
    "common_event_authorized",
    "common_event_completed",
    "diagnostic_fine_endpoint_is_accepted_state",
    "external_publication_authorized",
    "fine_path_commit_authorized",
    "fine_path_committed",
    "independent_method_agreement_earned",
    "journal",
    "ledger",
    "mechanism_result_authorized",
    "mechanism_result_earned",
    "new_protocol_id_authorized",
    "physical_result_authorized",
    "physical_result_earned",
    "physical_transition_claim_authorized",
    "production_SSPRK3_comparator_earned",
    "production_method_earned",
    "proto15_cursor",
    "push_authorized",
    "retained_EFT_evolution_authorized",
    "retry_4_authorized",
    "retry_5_authorized",
    "state_advance_authorized",
    "successor_remedy_selected",
    "temporal_retry_admission_called",
)
_RATIONAL_KEYS = frozenset({"numerator", "denominator"})
_DIFFERENCE_KEYS = frozenset(
    {
        "candidate_count",
        "co_maximizer_count",
        "coefficient_stream_sha256",
        "independent_evaluator_id",
        "independent_stationary_count_stream_sha256",
        "localization_classification",
        "lower",
        "maximum_candidates",
        "polynomial_count",
        "primary_evaluator_id",
        "primary_stationary_count_stream_sha256",
        "routes_agree",
        "survivor_key_stream_sha256",
        "upper",
    }
)
_CHANNEL_REPLAY_KEYS = frozenset(
    {
        "D01",
        "D12",
        "admission_passed",
        "candidate_evidence_available",
        "channel",
        "classification",
        "combined_coefficient_stream_sha256",
        "evaluator_id",
        "hash_evidence_available",
        "maximum_candidates_D01",
        "maximum_candidates_D12",
        "order_threshold_passed",
        "order_threshold_resolved",
        "refinement_depth",
        "route_evidence_available",
        "row_count",
        "row_stream_sha256",
        "sufficient_condition",
        "sufficient_condition_holds",
        "sufficient_contraction_failure",
        "sufficient_contraction_pass",
        "sufficient_pass_left",
        "sufficient_pass_right",
        "temporal_retry_permitted",
        "threshold_inconclusive",
    }
)
_COMPLETED_EXTRA_KEYS = frozenset(
    {
        "SSPRK3_stage_and_endpoint_record_count",
        "diagnostic_fine_endpoint",
        "exact_complete_C",
        "replay_receipt",
        "shadow_path_count",
        "shadow_proposal_count",
    }
)
_LOCALIZATION_CLASSES = frozenset(
    {"unique_maximum", "nonunique_or_interval_inconclusive"}
)
_MAX_RAW_LEAF_BYTES = 8 * 1024 * 1024
_MAX_STORE_LEAF_BYTES = 128 * 1024 * 1024
_MAX_PAYLOAD_MEMBER_BYTES = 32 * 1024 * 1024
_HEX = frozenset("0123456789abcdef")
_AUTHORITY_BLOBS = (
    (
        "Makefile",
        "ee9ffbe0edef2a9365ce0d3fadb8c6d7e4f6bcff2b552746d863479f3d389b4e",
    ),
    (
        "README.md",
        "33080916f125a1e8042bf40abef1c3ee4dc9066ee4dc559b6230f1becff10ab4",
    ),
    (
        "configs/fgc/fgc-1-tdg10-qa1-frz1.toml",
        "0d171c3caea98c4efa9f44c90e666d8a665ecc99ab273f2470d57e5d297f64b3",
    ),
    (
        "docs/claim-ledger.md",
        "d43d03a3899eaef883633c0552f5e6c7ab11429a9fc35005cbd1c9c7739d7174",
    ),
    (
        "docs/fgc-runtime-matrix.md",
        "3b65c09b916f9ce0751b4ce731d7b575ac38c073fa30622dd361b314660205c0",
    ),
    (
        "docs/fgc-tdg10-qa1-frz1.md",
        "d27e5a1f3b3f22ab017be0df3aa34a68410151f82977a8d8664a94650ac55b14",
    ),
    (
        "docs/research-roadmap.md",
        "da59eb2f0f8e512b2ed06e45f30fe8a200731913b527dd7700b5ad10240acbc0",
    ),
    (
        "results/README.md",
        "2dbdf191a16085d964162d85ae91d1ba68bd0c429cab661f1ce52ac1b181139a",
    ),
    (
        "results/fgc-1-tdg10-qa1-frz1.json",
        "87fd4fd750c176e17dbcdf2abeac48a7faa4d3e28f4071ba191b85103f9463ef",
    ),
    (
        "scripts/check_repo.py",
        "e5d53694bce58c8624241cd947fe6c4f662a535242c92b95660c7110385a48a5",
    ),
    (
        "scripts/reproduce_fgc_tdg10_qa1_frz1.py",
        "8d7849e3a80ad0041ece1228f727316089bf123a51bb8e9de668a6bd1ae6f876",
    ),
    (
        "scripts/run_fgc_tdg10_qa1.py",
        "036097ee42fba5602df521bac8b3fbd50100c9beeecd9ee118bdcfbe390a0576",
    ),
    (
        "src/recursive_horizons/fgc/evolution/tdg10_exact_complete_c_admission.py",
        "16d3f231c90eaa0b769a680cef0ae555b0f29fb3e2b9c551c9b51f3a290358af",
    ),
    (
        "src/recursive_horizons/fgc/evolution/tdg10_exact_complete_c_runtime.py",
        "8b4be91addd5d5927a63bc55e759f77369e3a29dabc8f0d350b73ab2b39ee352",
    ),
    (
        "src/recursive_horizons/fgc/evolution/tdg10_qa1_authority.py",
        "c2b164b449afb427a3d162aa2284d9b39777d8fdd5ff9b99cf5bd79a249a7b4e",
    ),
    (
        "tests/test_check_repo_tdg10_qa1_frz1.py",
        "ffb8e7ed7c4bd27c8cc8b6078229673317604676b38aa2c2097696848dbcfed8",
    ),
    (
        "tests/test_fgc_tdg10_exact_complete_c_admission.py",
        "2f0c7b006ec8d77c04ce3f61a633ed661b306a97926f416161d582bcf56bf946",
    ),
    (
        "tests/test_fgc_tdg10_exact_complete_c_runtime.py",
        "ac49d4c3925fa73bbd92a900d46b5159fe2265a39ee3ecbbe11ac635f3de10df",
    ),
    (
        "tests/test_fgc_tdg10_qa1_authority.py",
        "e326efb463f82262eaaa9061ec1e99ffe357d9e46510e68d13a6c579d079a2ee",
    ),
    (
        "tests/test_fgc_tdg10_qa1_runner.py",
        "9934a5f9b6770d84a7c1f3e5be1cf684a1ab3173ec5628723f3de8ee9f0e0a30",
    ),
)


class TDG10QA1PREF1Error(ValueError):
    """The compact contract, raw terminal, store, or independent replay differs."""

    def __init__(self, stop_id: str, detail: object) -> None:
        self.stop_id = str(stop_id)
        self.detail = " ".join(str(detail).split())[:640]
        super().__init__(f"{self.stop_id}: {self.detail}")


def _stop(stop_id: str, detail: object) -> NoReturn:
    raise TDG10QA1PREF1Error(stop_id, detail)


def _canonical(value: object) -> bytes:
    try:
        return json.dumps(
            value,
            sort_keys=True,
            separators=(",", ":"),
            ensure_ascii=True,
            allow_nan=False,
        ).encode("ascii")
    except (TypeError, ValueError) as exc:
        raise TDG10QA1PREF1Error("PREF1_CANONICAL_DRIFT", exc) from exc


def canonical_result(value: object) -> bytes:
    try:
        return (
            json.dumps(
                value,
                sort_keys=True,
                indent=2,
                ensure_ascii=True,
                allow_nan=False,
            )
            + "\n"
        ).encode("ascii")
    except (TypeError, ValueError) as exc:
        raise TDG10QA1PREF1Error("PREF1_COMPACT_DRIFT", exc) from exc


def _pairs(items: list[tuple[str, Any]]) -> dict[str, Any]:
    answer: dict[str, Any] = {}
    for key, value in items:
        if key in answer:
            raise ValueError(f"duplicate key {key}")
        answer[key] = value
    return answer


def _json(raw: bytes, label: str) -> dict[str, Any]:
    try:
        value = json.loads(
            raw.decode("ascii"),
            object_pairs_hook=_pairs,
            parse_constant=lambda item: (_ for _ in ()).throw(ValueError(item)),
        )
    except (UnicodeDecodeError, json.JSONDecodeError, ValueError) as exc:
        raise TDG10QA1PREF1Error("PREF1_JSON_DRIFT", label) from exc
    if not isinstance(value, dict) or canonical_result(value) != raw:
        _stop("PREF1_JSON_DRIFT", f"{label} is not canonical pretty JSON")
    return value


def _git(root: Path, *arguments: str) -> bytes:
    environment = {
        key: value for key, value in os.environ.items() if not key.startswith("GIT_")
    }
    environment.update(
        {
            "GIT_OPTIONAL_LOCKS": "0",
            "LC_ALL": "C",
            "LANG": "C",
        }
    )
    try:
        return subprocess.run(
            (
                "git",
                "--no-replace-objects",
                "--no-optional-locks",
                "-c",
                "core.fsmonitor=false",
                "-c",
                "core.untrackedCache=false",
                *arguments,
            ),
            cwd=root,
            env=environment,
            check=True,
            stdin=subprocess.DEVNULL,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
        ).stdout
    except (OSError, subprocess.CalledProcessError) as exc:
        raise TDG10QA1PREF1Error("PREF1_GIT_DRIFT", arguments) from exc


def _identity(metadata: os.stat_result) -> tuple[int, int, int, int, int]:
    return (
        metadata.st_dev,
        metadata.st_ino,
        metadata.st_size,
        metadata.st_mtime_ns,
        metadata.st_ctime_ns,
    )


def _open_directory(root: Path, relative: str) -> int:
    candidate = Path(relative)
    if candidate.is_absolute() or any(
        part in {"", ".", ".."} for part in candidate.parts
    ):
        _stop("PREF1_PATH_UNSAFE", relative)
    root_before = root.lstat()
    if stat.S_ISLNK(root_before.st_mode) or not stat.S_ISDIR(root_before.st_mode):
        _stop("PREF1_PATH_UNSAFE", root)
    descriptor = os.open(
        root,
        os.O_RDONLY | getattr(os, "O_DIRECTORY", 0) | getattr(os, "O_NOFOLLOW", 0),
    )
    try:
        if _identity(os.fstat(descriptor)) != _identity(root_before):
            _stop("PREF1_PATH_RACED", root)
        for part in candidate.parts:
            before = os.stat(part, dir_fd=descriptor, follow_symlinks=False)
            if stat.S_ISLNK(before.st_mode) or not stat.S_ISDIR(before.st_mode):
                _stop("PREF1_PATH_UNSAFE", relative)
            child = os.open(
                part,
                os.O_RDONLY
                | getattr(os, "O_DIRECTORY", 0)
                | getattr(os, "O_NOFOLLOW", 0),
                dir_fd=descriptor,
            )
            active = os.fstat(child)
            after = os.stat(part, dir_fd=descriptor, follow_symlinks=False)
            if (
                not stat.S_ISDIR(active.st_mode)
                or _identity(before) != _identity(active)
                or _identity(before) != _identity(after)
            ):
                os.close(child)
                _stop("PREF1_PATH_RACED", relative)
            os.close(descriptor)
            descriptor = child
        return descriptor
    except Exception:
        os.close(descriptor)
        raise


def _read_leaf_at(directory_fd: int, name: str, *, maximum: int, label: str) -> bytes:
    if not name or "/" in name or name in {".", ".."}:
        _stop("PREF1_PATH_UNSAFE", label)
    before = os.stat(name, dir_fd=directory_fd, follow_symlinks=False)
    if (
        stat.S_ISLNK(before.st_mode)
        or not stat.S_ISREG(before.st_mode)
        or before.st_nlink != 1
        or before.st_size < 0
        or before.st_size > maximum
    ):
        _stop("PREF1_LEAF_UNSAFE", label)
    descriptor = os.open(
        name,
        os.O_RDONLY | getattr(os, "O_NOFOLLOW", 0),
        dir_fd=directory_fd,
    )
    try:
        active = os.fstat(descriptor)
        if _identity(active) != _identity(before):
            _stop("PREF1_LEAF_RACED", label)
        chunks: list[bytes] = []
        remaining = before.st_size
        while remaining:
            block = os.read(descriptor, min(1 << 20, remaining))
            if not block:
                _stop("PREF1_LEAF_SHORT_READ", label)
            chunks.append(block)
            remaining -= len(block)
        if os.read(descriptor, 1):
            _stop("PREF1_LEAF_GREW", label)
        after = os.fstat(descriptor)
    finally:
        os.close(descriptor)
    final = os.stat(name, dir_fd=directory_fd, follow_symlinks=False)
    if _identity(before) != _identity(after) or _identity(before) != _identity(final):
        _stop("PREF1_LEAF_CHANGED", label)
    return b"".join(chunks)


def _read_repo_leaf(root: Path, relative: str, maximum: int) -> bytes:
    candidate = Path(relative)
    parent = candidate.parent.as_posix()
    descriptor = (
        os.open(
            root,
            os.O_RDONLY | getattr(os, "O_DIRECTORY", 0) | getattr(os, "O_NOFOLLOW", 0),
        )
        if parent == "."
        else _open_directory(root, parent)
    )
    try:
        return _read_leaf_at(
            descriptor, candidate.name, maximum=maximum, label=relative
        )
    finally:
        os.close(descriptor)


def _read_named_child(
    directory_fd: int,
    directory_name: str,
    expected_name: str,
    *,
    relative: str,
    maximum: int,
) -> bytes:
    edge_before = os.stat(directory_name, dir_fd=directory_fd, follow_symlinks=False)
    if stat.S_ISLNK(edge_before.st_mode) or not stat.S_ISDIR(edge_before.st_mode):
        _stop("PREF1_RAW_TREE_DRIFT", relative)
    child = os.open(
        directory_name,
        os.O_RDONLY | getattr(os, "O_DIRECTORY", 0) | getattr(os, "O_NOFOLLOW", 0),
        dir_fd=directory_fd,
    )
    try:
        before = os.fstat(child)
        if _identity(edge_before) != _identity(before):
            _stop("PREF1_RAW_TREE_CHANGED", relative)
        with os.scandir(child) as entries:
            names = tuple(sorted(item.name for item in entries))
        if names != (expected_name,):
            _stop("PREF1_RAW_TREE_DRIFT", (relative, names))
        raw = _read_leaf_at(child, expected_name, maximum=maximum, label=relative)
        with os.scandir(child) as entries:
            final_names = tuple(sorted(item.name for item in entries))
        after = os.fstat(child)
    finally:
        os.close(child)
    edge_after = os.stat(directory_name, dir_fd=directory_fd, follow_symlinks=False)
    if (
        names != final_names
        or _identity(before) != _identity(after)
        or _identity(edge_before) != _identity(edge_after)
    ):
        _stop("PREF1_RAW_TREE_CHANGED", relative)
    return raw


def _raw_tree(root: Path) -> dict[str, bytes]:
    descriptor = _open_directory(root, RAW_NAMESPACE)
    try:
        before = os.fstat(descriptor)
        with os.scandir(descriptor) as entries:
            names = tuple(sorted(item.name for item in entries))
        expected_names = (
            "fine-endpoint",
            "manifest.json",
            "payloads",
            "terminal.json",
        )
        if names != expected_names:
            _stop("PREF1_RAW_TREE_DRIFT", names)
        blobs = {
            "manifest.json": _read_leaf_at(
                descriptor,
                "manifest.json",
                maximum=1 << 20,
                label=f"{RAW_NAMESPACE}/manifest.json",
            ),
            "terminal.json": _read_leaf_at(
                descriptor,
                "terminal.json",
                maximum=_MAX_RAW_LEAF_BYTES,
                label=f"{RAW_NAMESPACE}/terminal.json",
            ),
            DESCRIPTOR_RELATIVE: _read_named_child(
                descriptor,
                "fine-endpoint",
                f"{DESCRIPTOR_CONTENT_ADDRESS}.json",
                relative=f"{RAW_NAMESPACE}/{DESCRIPTOR_RELATIVE}",
                maximum=_MAX_RAW_LEAF_BYTES,
            ),
            PAYLOAD_RELATIVE: _read_named_child(
                descriptor,
                "payloads",
                f"{PAYLOAD_SEMANTIC_SHA256}.npz",
                relative=f"{RAW_NAMESPACE}/{PAYLOAD_RELATIVE}",
                maximum=_MAX_RAW_LEAF_BYTES,
            ),
        }
        with os.scandir(descriptor) as entries:
            final_names = tuple(sorted(item.name for item in entries))
        after = os.fstat(descriptor)
    finally:
        os.close(descriptor)
    if names != final_names or _identity(before) != _identity(after):
        _stop("PREF1_RAW_TREE_CHANGED", RAW_NAMESPACE)
    observed = {
        "manifest.json": sha256(blobs["manifest.json"]).hexdigest(),
        "terminal.json": sha256(blobs["terminal.json"]).hexdigest(),
        DESCRIPTOR_RELATIVE: sha256(blobs[DESCRIPTOR_RELATIVE]).hexdigest(),
        PAYLOAD_RELATIVE: sha256(blobs[PAYLOAD_RELATIVE]).hexdigest(),
    }
    expected = {
        "manifest.json": RAW_MANIFEST_SHA256,
        "terminal.json": RAW_TERMINAL_SHA256,
        DESCRIPTOR_RELATIVE: DESCRIPTOR_RAW_SHA256,
        PAYLOAD_RELATIVE: PAYLOAD_RAW_SHA256,
    }
    if observed != expected:
        _stop("PREF1_RAW_HASH_DRIFT", observed)
    return blobs


def _raw_snapshot(
    root: Path,
) -> tuple[dict[str, bytes], dict[str, Any], dict[str, Any], dict[str, Any]]:
    blobs = _raw_tree(root)
    return (
        blobs,
        _json(blobs["manifest.json"], "raw manifest"),
        _json(blobs["terminal.json"], "raw terminal"),
        _json(blobs[DESCRIPTOR_RELATIVE], "diagnostic descriptor"),
    )


def _snapshot_store(root: Path) -> tuple[int, str]:
    """Reproduce the sealed stack-order store digest with no-following reads."""

    base_fd = _open_directory(root, STORE_PATH)
    digest = sha256(b"TDG9-LOC1-STORE-SNAPSHOT-v1\n")
    count = 0
    stack: list[tuple[int, str]] = [(base_fd, "")]
    try:
        while stack:
            directory_fd, prefix = stack.pop()
            try:
                before = os.fstat(directory_fd)
                with os.scandir(directory_fd) as entries:
                    items = sorted(entries, key=lambda item: item.name)
                initial_names = tuple(item.name for item in items)
                children: list[tuple[int, str]] = []
                for item in items:
                    metadata = item.stat(follow_symlinks=False)
                    relative = f"{prefix}/{item.name}" if prefix else item.name
                    if item.is_symlink():
                        _stop("PREF1_STORE_SYMLINK", relative)
                    if stat.S_ISDIR(metadata.st_mode):
                        child_fd = os.open(
                            item.name,
                            os.O_RDONLY
                            | getattr(os, "O_DIRECTORY", 0)
                            | getattr(os, "O_NOFOLLOW", 0),
                            dir_fd=directory_fd,
                        )
                        active = os.fstat(child_fd)
                        if _identity(metadata) != _identity(active):
                            os.close(child_fd)
                            _stop("PREF1_STORE_DIRECTORY_RACED", relative)
                        children.append((child_fd, relative))
                    elif stat.S_ISREG(metadata.st_mode):
                        raw = _read_leaf_at(
                            directory_fd,
                            item.name,
                            maximum=_MAX_STORE_LEAF_BYTES,
                            label=relative,
                        )
                        digest.update(
                            (relative + "\0" + sha256(raw).hexdigest() + "\n").encode(
                                "ascii"
                            )
                        )
                        count += 1
                    else:
                        _stop("PREF1_STORE_FOREIGN_TYPE", relative)
                with os.scandir(directory_fd) as entries:
                    final_names = tuple(sorted(item.name for item in entries))
                after = os.fstat(directory_fd)
                if initial_names != final_names or _identity(before) != _identity(
                    after
                ):
                    for child_fd, _ in children:
                        os.close(child_fd)
                    _stop("PREF1_STORE_DIRECTORY_CHANGED", prefix or ".")
                stack.extend(children)
            finally:
                os.close(directory_fd)
    except Exception:
        for directory_fd, _ in stack:
            try:
                os.close(directory_fd)
            except OSError:
                pass
        raise
    return count, digest.hexdigest()


def _mapping(value: object, keys: frozenset[str], label: str) -> Mapping[str, Any]:
    if not isinstance(value, Mapping) or set(value) != keys:
        _stop("PREF1_SCHEMA_DRIFT", (label, sorted(keys)))
    return value


def _digest(value: object, *, label: str) -> str:
    if (
        not isinstance(value, str)
        or len(value) != 64
        or any(character not in _HEX for character in value)
    ):
        _stop("PREF1_DIGEST_DRIFT", label)
    return value


def _count(value: object, *, label: str, allow_zero: bool = False) -> int:
    if (
        type(value) is not int
        or isinstance(value, bool)
        or value < (0 if allow_zero else 1)
    ):
        _stop("PREF1_COUNT_DRIFT", label)
    return value


def _encode_rational(value: Fraction) -> dict[str, str]:
    if not isinstance(value, Fraction) or isinstance(value, bool):
        _stop("PREF1_EXACT_TYPE_DRIFT", type(value).__name__)
    encoded = {
        "denominator": str(value.denominator),
        "numerator": str(value.numerator),
    }
    rebuilt = Fraction(int(encoded["numerator"]), int(encoded["denominator"]))
    if rebuilt != value or str(rebuilt.numerator) != encoded["numerator"]:
        _stop("PREF1_EXACT_ROUNDTRIP_DRIFT", encoded)
    return encoded


def _fraction_from_rational(value: object, *, label: str) -> Fraction:
    encoded = _mapping(value, _RATIONAL_KEYS, label)
    numerator = encoded["numerator"]
    denominator = encoded["denominator"]
    if not isinstance(numerator, str) or not isinstance(denominator, str):
        _stop("PREF1_EXACT_TYPE_DRIFT", label)
    try:
        fraction = Fraction(int(numerator), int(denominator))
    except (ValueError, ZeroDivisionError) as exc:
        raise TDG10QA1PREF1Error("PREF1_EXACT_PARSE_DRIFT", label) from exc
    if str(fraction.numerator) != numerator or str(fraction.denominator) != denominator:
        _stop("PREF1_EXACT_ROUNDTRIP_DRIFT", label)
    return fraction


def reduce_qa1_terminal(
    *,
    complete_admission_passed: bool,
    failed_channels: Sequence[str],
) -> str:
    if complete_admission_passed is True and tuple(failed_channels) == ():
        return RAW_CLASSIFICATION
    return NONPASS_CLASS


def _sealed_store() -> dict[str, Any]:
    return {
        "leaf_count": SEALED_STORE_LEAF_COUNT,
        "sha256": SEALED_STORE_SHA256,
    }


def _nonclaims() -> dict[str, Any]:
    return {
        "FGCQR_authorized": False,
        "GR0_calibration_authorized": False,
        "GR0_calibration_completed": False,
        "PDE_state_commit_authorized": False,
        "PDE_state_committed": False,
        "SGBL_authorized": False,
        "campaign_state_write_authorized": False,
        "candidate_branch_opened": False,
        "candidate_execution_authorized": False,
        "common_event_authorized": False,
        "common_event_completed": False,
        "diagnostic_fine_endpoint_is_accepted_state": False,
        "diagnostic_qualification_only": True,
        "external_publication_authorized": False,
        "fine_path_commit_authorized": False,
        "fine_path_committed": False,
        "independent_method_agreement_earned": False,
        "mechanism_result_authorized": False,
        "mechanism_result_earned": False,
        "new_protocol_id_authorized": False,
        "physical_result_authorized": False,
        "physical_result_earned": False,
        "physical_transition_claim_authorized": False,
        "production_SSPRK3_comparator_earned": False,
        "production_method_earned": False,
        "push_authorized": False,
        "retained_EFT_evolution_authorized": False,
        "retry_4_authorized": False,
        "retry_5_authorized": False,
        "state_advance_authorized": False,
        "successor_remedy_selected": False,
        "target_protocol": TARGET_PROTOCOL,
        "temporal_retry_admission_called": False,
    }


def _expected_config() -> dict[str, Any]:
    return {
        "schema_version": SCHEMA_VERSION,
        "artifact_id": ARTIFACT_ID,
        "classification": CLASSIFICATION,
        "target_protocol": TARGET_PROTOCOL,
        "owner_document": OWNER_DOCUMENT,
        "raw": {
            "namespace": RAW_NAMESPACE,
            "schema": RAW_SCHEMA,
            "classification": RAW_CLASSIFICATION,
            "runner_id": RUNNER_ID,
            "authority_commit": AUTHORITY_COMMIT,
            "manifest_sha256": RAW_MANIFEST_SHA256,
            "terminal_sha256": RAW_TERMINAL_SHA256,
            "descriptor_raw_sha256": DESCRIPTOR_RAW_SHA256,
            "descriptor_content_address": DESCRIPTOR_CONTENT_ADDRESS,
            "payload_raw_sha256": PAYLOAD_RAW_SHA256,
            "payload_semantic_sha256": PAYLOAD_SEMANTIC_SHA256,
            "leaf_count": 4,
        },
        "authority": {
            "parent_commit": AUTHORITY_PARENT,
            "blob_count": len(_AUTHORITY_BLOBS),
            "blobs": [
                {"path": path, "sha256": digest} for path, digest in _AUTHORITY_BLOBS
            ],
        },
        "store": {
            "path": STORE_PATH,
            "leaf_count": SEALED_STORE_LEAF_COUNT,
            "snapshot_sha256": SEALED_STORE_SHA256,
            "snapshot_algorithm": "TDG9-LOC1-STORE-SNAPSHOT-v1",
        },
        "selection": {
            "experiment_label": EXPERIMENT_LABEL,
            "tableau_selector": "SSPRK3",
            "tableau_runtime_selector": TABLEAU_RUNTIME_SELECTOR,
            "actual_spatial_operator": "inherited_RK4_2049_SBP4",
            "interval_owner": INTERVAL_OWNER,
            "classifier_owner": "classify_tdg6_channel",
            "production_SSPRK3_comparator": False,
            "independent_method_agreement": False,
            "member_key": MEMBER_KEY,
            "physical_state_sha256": PHYSICAL_STATE_SHA256,
            "accepted_time_hex": ACCEPTED_TIME_HEX,
            "point_count": POINT_COUNT,
            "owned_row_count": OWNED_ROW_COUNT,
            "retry": RETRY,
            "prior_retry_count": PRIOR_RETRY_COUNT,
            "attempted_width_hex": WIDTH_HEX,
            "predecessor_generation": PREDECESSOR_GENERATION,
            "forbidden_retry3_generation": FORBIDDEN_RETRY3_GENERATION,
            "checkpoint_journal_sequence": CHECKPOINT_JOURNAL_SEQUENCE,
            "historical_retry3_rejection_sequence": (
                HISTORICAL_RETRY3_REJECTION_SEQUENCE
            ),
            "per_channel_maximum_candidates_D01": (PER_CHANNEL_MAXIMUM_CANDIDATES_D01),
            "per_channel_maximum_candidates_D12": (PER_CHANNEL_MAXIMUM_CANDIDATES_D12),
            "primary_refinement_depth": PRIMARY_REFINEMENT_DEPTH,
            "channel_order": list(CHANNEL_ORDER),
            "UR1_PREF1_result_sha256": UR1_PREF1_RESULT_SHA256,
        },
        "work_budget": dict(WORK_BUDGET),
        "decision": {
            "complete_admission_is_all_of": True,
            "interval_owner": INTERVAL_OWNER,
            "classifier_owner": "classify_tdg6_channel",
            "sufficient_condition": SUFFICIENT_CONDITION,
            "equality_passes": True,
            "retry3_restores_generation": PREDECESSOR_GENERATION,
            "retry3_never_restores_generation_10": True,
            "checkpoint_journal_is_sequence_10": True,
            "historical_retry3_rejection_is_sequence_11": True,
            "diagnostic_fine_endpoint_is_accepted_state": False,
            "later_state_adoption_requires_separate_freeze": True,
            "successor_remedy_selected": False,
        },
        "scope": {
            "one_time_live_raw_and_store_authentication": True,
            "ordinary_verifier_raw_blind": True,
            "ordinary_verifier_store_blind": True,
            "ordinary_verifier_shadow_blind": True,
            "QA1_runner_decision_code_imported": False,
            "TI2_runner_imported": False,
            "QA1_authority_imported": False,
            "QA1_exact_runtime_imported": False,
            "QA1_exact_admission_imported": False,
            "lower_level_exact_arithmetic_and_dual_localizers_reused": True,
            "generation_10_restored": False,
            "seven_proposals_reconstructed": True,
            "production_SSPRK3_comparator": False,
            "independent_method_agreement": False,
            "campaign_store_mutation_authorized": False,
            "raw_namespace_mutation_authorized": False,
            "temporal_retry_admission_authorized": False,
            "PDE_state_commit_authorized": False,
            "successor_remedy_selected": False,
            "common_event_authorized": False,
            "GR0_calibration_authorized": False,
            "candidate_branches_authorized": False,
            "mechanism_result_authorized": False,
            "physical_result_authorized": False,
            "SGBL_authorized": False,
            "FGCQR_authorized": False,
            "retry_4_authorized": False,
            "retry_5_authorized": False,
        },
        "claims": {
            "QA1_raw_result_independently_bound": True,
            "all_eighteen_channels_admitted": True,
            "exact_complete_C_all_channel_pass": True,
            "diagnostic_fine_endpoint_only": True,
            "accepted_state": False,
            "production_method_earned": False,
            "production_SSPRK3_comparator": False,
            "independent_method_agreement": False,
            "state_advance_authorized": False,
            "successor_remedy_selected": False,
            "common_event_completed": False,
            "GR0_calibration_completed": False,
            "candidate_execution_authorized": False,
            "mechanism_result_earned": False,
            "physical_result_earned": False,
            "retained_EFT_evolution_authorized": False,
            "physical_transition_claim_authorized": False,
            "SGBL_authorized": False,
            "FGCQR_authorized": False,
            "external_publication_authorized": False,
        },
    }


def expected_config() -> dict[str, Any]:
    """Return the exact typed config contract this binder will accept.

    The checked-in scaffold was intentionally allowed to arrive before the
    binder.  Exposing a detached copy makes reconciliation mechanical while
    keeping this module's twenty-blob authority and anti-common-mode scope the
    source of truth.
    """

    return json.loads(_canonical(_expected_config()).decode("ascii"))


def _config(raw: bytes) -> dict[str, Any]:
    try:
        value = tomllib.loads(raw.decode("utf-8"))
    except (UnicodeDecodeError, tomllib.TOMLDecodeError) as exc:
        raise TDG10QA1PREF1Error("PREF1_CONFIG_DRIFT", exc) from exc
    if value != _expected_config():
        _stop("PREF1_CONFIG_DRIFT", "typed PREF1 config differs")
    return value


def expected_evidence(config_raw: bytes) -> dict[str, Any]:
    config = _config(config_raw)
    return {
        "schema_version": SCHEMA_VERSION,
        "artifact_id": ARTIFACT_ID,
        "classification": CLASSIFICATION,
        "target_protocol": TARGET_PROTOCOL,
        "raw": dict(config["raw"]),
        "authority": dict(config["authority"]),
        "store": dict(config["store"]),
        "selection": dict(config["selection"]),
        "work_budget": dict(config["work_budget"]),
        "decision": dict(config["decision"]),
        "scope": dict(config["scope"]),
        "claims": dict(config["claims"]),
    }


def _expected_channel_summaries() -> list[dict[str, Any]]:
    return [
        {
            "admission_passed": True,
            "channel": channel,
            "classification": "resolved_order_pass",
            "D12_upper": {
                "denominator": denominator,
                "numerator": numerator,
            },
            "order_threshold_passed": True,
            "order_threshold_resolved": True,
            "sufficient_condition": SUFFICIENT_CONDITION,
            "sufficient_contraction_pass": True,
            "temporal_retry_permitted": False,
        }
        for channel, numerator, denominator in EXPECTED_D12_UPPERS
    ]


def _expected_independent_replay_header() -> dict[str, Any]:
    return {
        "SSPRK3_stage_and_endpoint_record_count": 28,
        "admitted_count": 18,
        "all_of_identity_verified": True,
        "channel_count": 18,
        "complete_admission_passed": True,
        "failed_channels": [],
        "failed_count": 0,
        "independent_route_required": True,
        "interval_owner": INTERVAL_OWNER,
        "predecessor_generation": PREDECESSOR_GENERATION,
        "raw_terminal_class_matches_independent_reduction": True,
        "replay_receipt": dict(EXPECTED_REPLAY_RECEIPT),
        "retry3_never_restores_generation_10": True,
        "shadow_path_count": 7,
        "shadow_proposal_count": 7,
        "sufficient_condition": SUFFICIENT_CONDITION,
    }


def _expected_diagnostic_binding() -> dict[str, Any]:
    return {
        "accepted_state": False,
        "content_address_verified": True,
        "descriptor_content_address": DESCRIPTOR_CONTENT_ADDRESS,
        "descriptor_raw_sha256": DESCRIPTOR_RAW_SHA256,
        "descriptor_relative": DESCRIPTOR_RELATIVE,
        "diagnostic_qualification_only": True,
        "every_state_calibration_candidate_physics_flag_false": True,
        "payload_manifest_verified": True,
        "payload_raw_sha256": PAYLOAD_RAW_SHA256,
        "payload_relative": PAYLOAD_RELATIVE,
        "payload_semantic_sha256": PAYLOAD_SEMANTIC_SHA256,
        "present": True,
    }


def _expected_conclusion() -> dict[str, Any]:
    return {
        "PDE_or_candidate_state_opened": False,
        "all_eighteen_exact_complete_C_channels_pass": True,
        "diagnostic_fine_endpoint_only": True,
        "independent_method_agreement": False,
        "licenses_only_later_separately_frozen_state_adoption": True,
        "physics_inference_permitted": False,
        "production_SSPRK3_comparator": False,
        "production_method_earned": False,
        "retry_4_authorized": False,
        "retry_5_authorized": False,
        "SGBL_authorized": False,
        "FGCQR_authorized": False,
        "state_advance_authorized": False,
        "successor_remedy_selected": False,
        "accepted_state": False,
        "common_event_completed": False,
        "GR0_calibration_completed": False,
        "external_publication_authorized": False,
    }


def expected_bound_payload(config_raw: bytes) -> dict[str, Any]:
    expected = expected_evidence(config_raw)
    sealed = _sealed_store()
    return {
        **expected,
        "raw_binding": {
            "canonical_duplicate_free_schema_verified": True,
            "descriptor_content_address": DESCRIPTOR_CONTENT_ADDRESS,
            "descriptor_raw_sha256": DESCRIPTOR_RAW_SHA256,
            "leaf_count": 4,
            "manifest_sha256": RAW_MANIFEST_SHA256,
            "nested_four_leaf_tree_verified": True,
            "payload_raw_sha256": PAYLOAD_RAW_SHA256,
            "payload_semantic_sha256": PAYLOAD_SEMANTIC_SHA256,
            "raw_namespace_unchanged": True,
            "terminal_sha256": RAW_TERMINAL_SHA256,
        },
        "authority_binding": {
            "all_committed_blob_hashes_verified": True,
            "blob_count": len(_AUTHORITY_BLOBS),
            "commit": AUTHORITY_COMMIT,
            "parent_commit": AUTHORITY_PARENT,
        },
        "store_binding": {
            "binder_snapshot_after": dict(sealed),
            "binder_snapshot_before": dict(sealed),
            "raw_snapshot_after": dict(sealed),
            "raw_snapshot_before": dict(sealed),
            "store_unchanged": True,
        },
        "diagnostic_binding": _expected_diagnostic_binding(),
        "independent_replay": {
            **_expected_independent_replay_header(),
            "channel_summaries": _expected_channel_summaries(),
        },
        "conclusion": _expected_conclusion(),
    }


def expected_compact_result(config_raw: bytes) -> dict[str, Any]:
    return {
        "artifact_id": ARTIFACT_ID,
        "artifact_payload": expected_bound_payload(config_raw),
    }


def _authenticate_authority(root: Path, manifest: Mapping[str, Any]) -> None:
    if manifest.get("authority_commit") != AUTHORITY_COMMIT:
        _stop("PREF1_AUTHORITY_DRIFT", manifest.get("authority_commit"))
    try:
        top = Path(_git(root, "rev-parse", "--show-toplevel").decode("utf-8").strip())
        if not os.path.samefile(root, top):
            _stop("PREF1_AUTHORITY_DRIFT", "repository is not the worktree root")
    except (OSError, UnicodeDecodeError) as exc:
        raise TDG10QA1PREF1Error(
            "PREF1_AUTHORITY_DRIFT", "Git worktree root differs"
        ) from exc
    if _git(root, "cat-file", "-t", AUTHORITY_COMMIT).strip() != b"commit":
        _stop("PREF1_AUTHORITY_DRIFT", "authority object is not a commit")
    parents = _git(root, "show", "-s", "--format=%P", AUTHORITY_COMMIT).split()
    if parents != [AUTHORITY_PARENT.encode("ascii")]:
        _stop("PREF1_AUTHORITY_LINEAGE_DRIFT", parents)
    if len(_AUTHORITY_BLOBS) != 20:
        _stop("PREF1_AUTHORITY_BLOB_COUNT_DRIFT", len(_AUTHORITY_BLOBS))
    expected_paths = tuple(path for path, _digest_value in _AUTHORITY_BLOBS)
    if len(set(expected_paths)) != len(expected_paths):
        _stop("PREF1_AUTHORITY_BLOB_COUNT_DRIFT", "duplicate path")
    observed_paths = tuple(
        sorted(
            _git(
                root,
                "diff-tree",
                "--root",
                "--no-commit-id",
                "--name-only",
                "-r",
                AUTHORITY_COMMIT,
            )
            .decode("utf-8")
            .splitlines()
        )
    )
    if observed_paths != tuple(sorted(expected_paths)):
        _stop("PREF1_AUTHORITY_DELTA_DRIFT", observed_paths)
    for relative, expected in _AUTHORITY_BLOBS:
        observed = sha256(
            _git(root, "show", f"{AUTHORITY_COMMIT}:{relative}")
        ).hexdigest()
        if observed != expected:
            _stop("PREF1_AUTHORITY_BLOB_DRIFT", relative)


def _expected_manifest() -> dict[str, Any]:
    return {
        "actual_spatial_operator": "inherited_RK4_2049_SBP4",
        "artifact_id": QA1_ARTIFACT_ID,
        "authority_commit": AUTHORITY_COMMIT,
        "diagnostic_fine_endpoint_is_accepted_state": False,
        "diagnostic_qualification_only": True,
        "output_leaves": list(EXPECTED_OUTPUT_LEAVES),
        "predecessor_generation": PREDECESSOR_GENERATION,
        "retry": RETRY,
        "runner_id": RUNNER_ID,
        "schema": RAW_SCHEMA,
        "state_advance_authorized": False,
        "tableau_selector": "SSPRK3",
        "target_protocol": TARGET_PROTOCOL,
    }


def _expected_terminal_base() -> dict[str, Any]:
    sealed = _sealed_store()
    return {
        **_nonclaims(),
        "actual_spatial_operator": "inherited_RK4_2049_SBP4",
        "artifact_id": QA1_ARTIFACT_ID,
        "attempted_width_hex": WIDTH_HEX,
        "authority_commit": AUTHORITY_COMMIT,
        "channel_order": list(CHANNEL_ORDER),
        "checkpoint_journal_sequence": CHECKPOINT_JOURNAL_SEQUENCE,
        "classification": RAW_CLASSIFICATION,
        "historical_retry3_rejection_sequence": (HISTORICAL_RETRY3_REJECTION_SEQUENCE),
        "predecessor_generation": PREDECESSOR_GENERATION,
        "retry": RETRY,
        "runner_id": RUNNER_ID,
        "schema": RAW_SCHEMA,
        "store_snapshot_after": sealed,
        "store_snapshot_before": sealed,
        "store_unchanged": True,
        "tableau_runtime_selector": TABLEAU_RUNTIME_SELECTOR,
        "tableau_selector": "SSPRK3",
    }


def _class_name(value: object) -> str:
    kind = type(value)
    return f"{kind.__module__}.{kind.__qualname__}"


def _callable_name(value: object) -> str:
    return (
        f"{getattr(value, '__module__', type(value).__module__)}."
        f"{getattr(value, '__qualname__', type(value).__qualname__)}"
    )


def _fingerprint_exact(value: object) -> object:
    if isinstance(value, Fraction):
        return {
            "numerator": str(value.numerator),
            "denominator": str(value.denominator),
        }
    if type(value) is float:
        if not math.isfinite(value):
            _stop("PREF1_FINGERPRINT_DRIFT", "nonfinite binary64")
        return {"binary64_hex": value.hex()}
    if is_dataclass(value) and not isinstance(value, type):
        return _fingerprint_exact(asdict(value))
    if isinstance(value, Mapping):
        return {str(key): _fingerprint_exact(item) for key, item in value.items()}
    if isinstance(value, (tuple, list)):
        return [_fingerprint_exact(item) for item in value]
    if value is None or type(value) in {str, int, bool}:
        return value
    _stop("PREF1_FINGERPRINT_DRIFT", type(value).__name__)


def _tracer_history_hash(member: object) -> str:
    arrays = [
        np.asarray(member.tracers.labels, dtype=np.float64),
        np.asarray(member.tracers.positions, dtype=np.float64),
        np.asarray(member.tracers.proper_times, dtype=np.float64),
    ]
    arrays.extend(
        np.asarray(row, dtype=np.float64) for row in member.tracers.event_proper_times
    )
    arrays.extend(
        np.asarray(row, dtype=np.float64) for row in member.tracers.event_fields
    )
    return array_content_sha256(*arrays)


def _member_fingerprint(member: object, descriptor_sha256: str) -> dict[str, object]:
    ledger = member.temporal_ledger
    if ledger is None:
        _stop("PREF1_REPLAY_DRIFT", "temporal ledger absent")
    return {
        "member_key": member.key,
        "method_label": member.method_label,
        "source_integrator": member.integrator_id,
        "point_count": member.point_count,
        "accepted_time_hex": float(member.time).hex(),
        "state_sha256": array_content_sha256(
            member.state.u, member.state.p, member.state.q
        ),
        "coordinates_sha256": array_content_sha256(member.initial.grid.coordinates),
        "grid_spacing_hex": float(member.initial.grid.spacing).hex(),
        "outer_radius_hex": float(member.initial.grid.maximum).hex(),
        "descriptor_sha256": descriptor_sha256,
        "operator_class": _class_name(member.operator),
        "projector_callable": _callable_name(member.projector),
        "transaction_class": _class_name(member.transaction),
        "tracer_class": _class_name(member.tracers),
        "ledger_class": _class_name(ledger),
        "transaction_sha256": sha256(
            _canonical(
                _fingerprint_exact(
                    {
                        "monitor": asdict(member.transaction.state),
                        "causal": asdict(member.transaction.causal_state),
                        "ledger": asdict(ledger),
                        "step_index": member.step_index,
                        "transaction_serial": member.transaction_serial,
                    }
                )
            )
        ).hexdigest(),
        "tracer_history_sha256": _tracer_history_hash(member),
    }


def _journal_evidence(root: Path) -> Mapping[str, Any]:
    relative = (
        f"{STORE_PATH}/journal/"
        f"{HISTORICAL_RETRY3_REJECTION_SEQUENCE:020d}-"
        f"{HISTORICAL_RETRY3_REJECTION_SHA256}.journal"
    )
    raw = _read_repo_leaf(root, relative, 8 << 20)
    if sha256(raw).hexdigest() != HISTORICAL_RETRY3_REJECTION_RAW_SHA256:
        _stop("PREF1_JOURNAL_HASH_DRIFT", relative)
    value = json.loads(raw.decode("ascii"), object_pairs_hook=_pairs)
    if (
        not isinstance(value, Mapping)
        or value.get("sequence") != HISTORICAL_RETRY3_REJECTION_SEQUENCE
        or value.get("record_sha256") != HISTORICAL_RETRY3_REJECTION_SHA256
        or value.get("kind") != "tdg6_rejection"
    ):
        _stop("PREF1_JOURNAL_EVIDENCE_DRIFT", "historical rejection")
    evidence = value.get("payload", {}).get("evidence")
    if not isinstance(evidence, Mapping):
        _stop("PREF1_JOURNAL_EVIDENCE_DRIFT", "evidence absent")
    if sha256(_canonical(evidence)).hexdigest() != HISTORICAL_JOURNAL_SHA256:
        _stop("PREF1_HISTORICAL_JOURNAL_DRIFT", HISTORICAL_JOURNAL_SHA256)
    return evidence


def _verify_generation9_checkpoint(root: Path) -> None:
    relative = (
        f"{STORE_PATH}/checkpoints/"
        f"{PREDECESSOR_GENERATION:020d}-{CHECKPOINT_SHA256}.json"
    )
    raw = _read_repo_leaf(root, relative, 8 << 20)
    if sha256(raw).hexdigest() != CHECKPOINT_RAW_SHA256:
        _stop("PREF1_CHECKPOINT_RAW_DRIFT", relative)
    value = json.loads(raw.decode("ascii"), object_pairs_hook=_pairs)
    if (
        not isinstance(value, Mapping)
        or value.get("generation") != PREDECESSOR_GENERATION
        or value.get("generation") == FORBIDDEN_RETRY3_GENERATION
        or value.get("checkpoint_sha256") != CHECKPOINT_SHA256
        or value.get("journal_sequence") != CHECKPOINT_JOURNAL_SEQUENCE
        or value.get("journal_tip_sha256") != CHECKPOINT_JOURNAL_SHA256
    ):
        _stop("PREF1_CHECKPOINT_DRIFT", "generation-9 / sequence-10 tip")
    journal_relative = (
        f"{STORE_PATH}/journal/"
        f"{CHECKPOINT_JOURNAL_SEQUENCE:020d}-{CHECKPOINT_JOURNAL_SHA256}.journal"
    )
    journal_raw = _read_repo_leaf(root, journal_relative, 8 << 20)
    if sha256(journal_raw).hexdigest() != CHECKPOINT_JOURNAL_RAW_SHA256:
        _stop("PREF1_CHECKPOINT_JOURNAL_DRIFT", journal_relative)
    journal = json.loads(journal_raw.decode("ascii"), object_pairs_hook=_pairs)
    if (
        not isinstance(journal, Mapping)
        or journal.get("sequence") != CHECKPOINT_JOURNAL_SEQUENCE
        or journal.get("record_sha256") != CHECKPOINT_JOURNAL_SHA256
        or journal.get("kind") == "tdg6_rejection"
    ):
        _stop("PREF1_CHECKPOINT_JOURNAL_DRIFT", "sequence 10 is not the tip")


def _restore_shadow(
    root: Path,
) -> tuple[tdg6.TDG6PreparedGR0Compositor, dict[str, Any]]:
    replay = dict(RETRY3_REPLAY)
    generation = int(replay["predecessor_generation"])
    if (
        generation != PREDECESSOR_GENERATION
        or generation == FORBIDDEN_RETRY3_GENERATION
    ):
        _stop("PREF1_GENERATION_DRIFT", generation)
    _verify_generation9_checkpoint(root)
    _journal_evidence(root)
    store = HLT16CampaignStore(root / STORE_PATH)
    shells = build_static_gr0_shells(root)
    checkpoint = store.authenticated_checkpoint_at_generation(generation)
    if (
        checkpoint.sha256 != CHECKPOINT_SHA256
        or checkpoint.generation != PREDECESSOR_GENERATION
        or checkpoint.generation == FORBIDDEN_RETRY3_GENERATION
        or checkpoint.journal_sequence != CHECKPOINT_JOURNAL_SEQUENCE
    ):
        _stop("PREF1_CHECKPOINT_DRIFT", generation)
    state = checkpoint.members.get(MEMBER_KEY)
    if state is None:
        _stop("PREF1_MEMBER_ABSENT", MEMBER_KEY)
    member = deepcopy(shells[MEMBER_KEY])
    try:
        campaign_runtime.restore_member_with_overlay(
            store, checkpoint, member, key=MEMBER_KEY
        )
        cursor = p15.Proto15Cursor(dict(state.cursor))
        cursor.validate()
    except Exception as exc:
        raise TDG10QA1PREF1Error("PREF1_RESTORE_DRIFT", RETRY) from exc
    ledger = member.temporal_ledger
    if (
        ledger is None
        or member.key != MEMBER_KEY
        or member.method_label != "RK4"
        or member.integrator_id != PRIMARY_METHOD
        or member.point_count != POINT_COUNT
        or float(member.time).hex() != ACCEPTED_TIME_HEX
        or ledger.current_macro_step_temporal_retry_count != PRIOR_RETRY_COUNT
        or cursor.mode != "RETRY_PENDING"
        or state.pending_owner != "temporal"
    ):
        _stop("PREF1_PREDECESSOR_DRIFT", RETRY)
    fingerprint = _member_fingerprint(member, state.descriptor_sha256)
    fixed = {
        "member_key": MEMBER_KEY,
        "method_label": "RK4",
        "source_integrator": PRIMARY_METHOD,
        "point_count": POINT_COUNT,
        "accepted_time_hex": ACCEPTED_TIME_HEX,
        "state_sha256": PHYSICAL_STATE_SHA256,
        "coordinates_sha256": COORDINATES_SHA256,
        "grid_spacing_hex": GRID_SPACING_HEX,
        "outer_radius_hex": OUTER_RADIUS_HEX,
        "descriptor_sha256": MEMBER_DESCRIPTOR_SHA256,
        "transaction_sha256": TRANSACTION_SHA256,
    }
    if any(fingerprint.get(name) != expected for name, expected in fixed.items()):
        _stop("PREF1_FINGERPRINT_DRIFT", RETRY)
    before = dict(fingerprint)
    try:
        prepared = tdg6.prepare_tdg6_gr0_compositor(
            method=COMPARATOR_METHOD,
            time=member.time,
            step_size=float.fromhex(WIDTH_HEX),
            state=member.state,
            rhs=member.operator,
            projector=member.projector,
            transaction=member.transaction,
            tracers=member.tracers,
            coordinates=member.initial.grid.coordinates,
            temporal_ledger=ledger,
            previous_step_index=member.step_index,
            previous_transaction_serial=member.transaction_serial,
        )
    except Exception as exc:
        raise TDG10QA1PREF1Error("PREF1_SHADOW_PREMISE_STOP", RETRY) from exc
    after = _member_fingerprint(member, state.descriptor_sha256)
    paths = (prepared.outer, prepared.medium, prepared.fine)
    proposals = tuple(attempt.proposal for path in paths for attempt in path.attempts)
    if (
        after != before
        or prepared.method != COMPARATOR_METHOD
        or prepared.initial_state_sha256 != PHYSICAL_STATE_SHA256
        or tuple(len(path.attempts) for path in paths) != (1, 2, 4)
        or len(proposals) != 7
        or any(proposal.method != COMPARATOR_METHOD for proposal in proposals)
        or any(len(proposal.stages) != 4 for proposal in proposals)
    ):
        _stop("PREF1_ONE_VARIABLE_IDENTITY_DRIFT", RETRY)
    receipt = {
        "member_key": MEMBER_KEY,
        "retry": RETRY,
        "predecessor_generation": PREDECESSOR_GENERATION,
        "checkpoint_journal_sequence": CHECKPOINT_JOURNAL_SEQUENCE,
        "historical_retry3_rejection_sequence": (HISTORICAL_RETRY3_REJECTION_SEQUENCE),
        "attempted_width_hex": WIDTH_HEX,
        "accepted_time_hex": fingerprint["accepted_time_hex"],
        "state_sha256": fingerprint["state_sha256"],
        "descriptor_sha256": fingerprint["descriptor_sha256"],
        "transaction_sha256": fingerprint["transaction_sha256"],
        "historical_journal_sha256": HISTORICAL_JOURNAL_SHA256,
    }
    if receipt != EXPECTED_REPLAY_RECEIPT:
        _stop("PREF1_REPLAY_RECEIPT_DRIFT", receipt)
    return prepared, dict(EXPECTED_REPLAY_RECEIPT)


@dataclass(frozen=True, slots=True)
class _IndependentDifference:
    lower: Fraction
    upper: Fraction
    polynomial_count: int
    candidate_count: int
    co_maximizer_count: int
    localization_classification: str
    coefficient_stream_sha256: str
    survivor_key_stream_sha256: str
    primary_stationary_count_stream_sha256: str
    independent_stationary_count_stream_sha256: str
    primary_evaluator_id: str
    independent_evaluator_id: str
    maximum_candidates: int
    routes_agree: bool = True


@dataclass(frozen=True, slots=True)
class _IndependentAdmission:
    d01: _IndependentDifference
    d12: _IndependentDifference
    decision: object
    row_count: int
    row_stream_sha256: str
    combined_coefficient_stream_sha256: str
    sufficient_pass_left: Fraction
    sufficient_pass_right: Fraction
    sufficient_contraction_pass: bool
    sufficient_contraction_failure: bool
    threshold_inconclusive: bool
    maximum_candidates_D01: int
    maximum_candidates_D12: int
    refinement_depth: int
    evaluator_id: str = TDG10_EVALUATOR_ID


@dataclass(frozen=True, slots=True)
class _IndependentChannel:
    channel: str
    evidence: _IndependentAdmission


@dataclass(frozen=True, slots=True)
class _IndependentAssessment:
    channel_evidence: tuple[_IndependentChannel, ...]
    complete_admission_passed: bool
    failed_channels: tuple[str, ...]


def _binary64(value: object, *, label: str) -> float:
    if type(value) is float:
        answer = value
    elif type(value) is np.float64:
        answer = float(value)
    else:
        _stop("PREF1_BINARY64_DRIFT", label)
    if not math.isfinite(answer):
        _stop("PREF1_BINARY64_DRIFT", f"{label} is nonfinite")
    return answer


def _same_binary64(left: object, right: object) -> bool:
    first = _binary64(left, label="left binary64")
    second = _binary64(right, label="right binary64")
    return np.float64(first).tobytes() == np.float64(second).tobytes()


def _state_sha256(state: object) -> str:
    try:
        return array_content_sha256(state.u, state.p, state.q)
    except Exception as exc:
        raise TDG10QA1PREF1Error("PREF1_STATE_DRIFT", type(state).__name__) from exc


def _rhs_sha256(rhs: object) -> str:
    try:
        return array_content_sha256(rhs.du, rhs.dp, rhs.dq)
    except Exception as exc:
        raise TDG10QA1PREF1Error("PREF1_RHS_DRIFT", type(rhs).__name__) from exc


def _owned_state(state: object) -> np.ndarray:
    arrays = tuple(np.asarray(getattr(state, name)) for name in ("u", "p", "q"))
    if any(
        value.dtype != np.dtype("<f8")
        or value.shape != (POINT_COUNT, 6)
        or not value.flags.c_contiguous
        or not np.isfinite(value).all()
        for value in arrays
    ):
        _stop("PREF1_OWNED_STATE_DRIFT", [value.shape for value in arrays])
    answer = np.stack(tuple(value[1:-4] for value in arrays), axis=0)
    if answer.shape != (3, OWNED_ROW_COUNT, 6):
        _stop("PREF1_OWNED_STATE_DRIFT", answer.shape)
    return answer


def _owned_rhs(rhs: object) -> np.ndarray:
    arrays = tuple(np.asarray(getattr(rhs, name)) for name in ("du", "dp", "dq"))
    if any(
        value.dtype != np.dtype("<f8")
        or value.shape != (POINT_COUNT, 6)
        or not value.flags.c_contiguous
        or not np.isfinite(value).all()
        for value in arrays
    ):
        _stop("PREF1_OWNED_RHS_DRIFT", [value.shape for value in arrays])
    answer = np.stack(tuple(value[1:-4] for value in arrays), axis=0)
    if answer.shape != (3, OWNED_ROW_COUNT, 6):
        _stop("PREF1_OWNED_RHS_DRIFT", answer.shape)
    return answer


def _endpoint_rhs(proposal: object) -> tuple[object, object]:
    stages = tuple(getattr(proposal, "stages", ()))
    if len(stages) != 4 or tuple(
        getattr(item, "stage_name", None) for item in stages
    ) != ("ssprk3_s0", "ssprk3_s1", "ssprk3_s2", "candidate_endpoint"):
        _stop("PREF1_STAGE_RECORD_DRIFT", len(stages))
    first, last = stages[0], stages[-1]
    if (
        not _same_binary64(first.time, proposal.initial_time)
        or not _same_binary64(last.time, proposal.final_time)
        or _state_sha256(first.state) != _state_sha256(proposal.initial_state)
        or _state_sha256(last.state) != _state_sha256(proposal.candidate_state)
    ):
        _stop("PREF1_ENDPOINT_RHS_DRIFT", proposal.method)
    return first.rhs, last.rhs


def _validate_prepared_shadow(prepared: object) -> float:
    if getattr(prepared, "method", None) != COMPARATOR_METHOD:
        _stop("PREF1_TABLEAU_DRIFT", getattr(prepared, "method", None))
    start = _binary64(prepared.initial_time, label="prepared.initial_time")
    final = _binary64(prepared.final_time, label="prepared.final_time")
    width = final - start
    if not math.isfinite(width) or width <= 0.0:
        _stop("PREF1_SHADOW_WIDTH_DRIFT", width)
    if prepared.initial_state_sha256 != PHYSICAL_STATE_SHA256:
        _stop("PREF1_INITIAL_STATE_DRIFT", prepared.initial_state_sha256)
    initial_rhs: list[str] = []
    for path, level, count in (
        (prepared.outer, "outer", 1),
        (prepared.medium, "medium", 2),
        (prepared.fine, "fine", 4),
    ):
        attempts = tuple(path.attempts)
        if (
            path.level != level
            or path.method != COMPARATOR_METHOD
            or len(attempts) != count
            or path.initial_state_sha256 != PHYSICAL_STATE_SHA256
            or not _same_binary64(path.initial_time, start)
            or not _same_binary64(path.final_time, final)
        ):
            _stop("PREF1_SHADOW_PATH_DRIFT", level)
        step = width / count
        boundaries = tuple(start + index * step for index in range(count + 1))
        if not _same_binary64(boundaries[-1], final):
            _stop("PREF1_SHADOW_BOUNDARY_DRIFT", level)
        previous_state: str | None = None
        previous_rhs: str | None = None
        for index, (attempt, left, right) in enumerate(
            zip(attempts, boundaries[:-1], boundaries[1:], strict=True)
        ):
            if attempt.accepted is None or attempt.retry is not None:
                _stop("PREF1_SHADOW_ATTEMPT_DRIFT", (level, index))
            proposal = attempt.proposal
            left_rhs, right_rhs = _endpoint_rhs(proposal)
            if (
                proposal.method != COMPARATOR_METHOD
                or not _same_binary64(proposal.initial_time, left)
                or not _same_binary64(proposal.final_time, right)
                or not _same_binary64(right - left, step)
            ):
                _stop("PREF1_SHADOW_BOUNDARY_DRIFT", (level, index))
            state_digest = _state_sha256(proposal.initial_state)
            if index == 0:
                if state_digest != PHYSICAL_STATE_SHA256:
                    _stop("PREF1_SHADOW_INITIAL_DRIFT", level)
                initial_rhs.append(_rhs_sha256(left_rhs))
            elif (
                state_digest != previous_state or _rhs_sha256(left_rhs) != previous_rhs
            ):
                _stop("PREF1_SHADOW_CONTIGUITY_DRIFT", (level, index))
            previous_state = _state_sha256(proposal.candidate_state)
            previous_rhs = _rhs_sha256(right_rhs)
            if _state_sha256(attempt.accepted.state) != previous_state:
                _stop("PREF1_SHADOW_ACCEPTED_DRIFT", (level, index))
    if len(set(initial_rhs)) != 1:
        _stop("PREF1_SHARED_INITIAL_RHS_DRIFT", initial_rhs)
    _owned_state(prepared.outer.attempts[0].proposal.initial_state)
    return width


def _segment(
    proposal: object,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, float]:
    left_rhs, right_rhs = _endpoint_rhs(proposal)
    width = _binary64(proposal.final_time, label="proposal.final_time") - _binary64(
        proposal.initial_time, label="proposal.initial_time"
    )
    if not math.isfinite(width) or width <= 0.0:
        _stop("PREF1_SEGMENT_WIDTH_DRIFT", width)
    return (
        _owned_state(proposal.initial_state),
        _owned_rhs(left_rhs),
        _owned_state(proposal.candidate_state),
        _owned_rhs(right_rhs),
        width,
    )


def _path_segments(
    path: object,
) -> tuple[tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, float], ...]:
    return tuple(_segment(attempt.proposal) for attempt in path.attempts)


def _channel_rows(
    prepared: object, channel: str, *, width: float
) -> tuple[tuple[object, object, object], ...]:
    outer = _path_segments(prepared.outer)
    medium = _path_segments(prepared.medium)
    fine = _path_segments(prepared.fine)
    if len(outer) != 1 or len(medium) != 2 or len(fine) != 4:
        _stop("PREF1_SHADOW_COUNT_DRIFT", channel)
    if (
        not _same_binary64(outer[0][4], width)
        or any(not _same_binary64(item[4], width / 2.0) for item in medium)
        or any(not _same_binary64(item[4], width / 4.0) for item in fine)
    ):
        _stop("PREF1_SEGMENT_WIDTH_DRIFT", channel)
    block_name, field_name = channel.split(":", 1)
    block = ("u", "p", "q").index(block_name)
    field = ("alpha", "v", "lambda", "R", "phi", "chi").index(field_name)

    def scalar(
        item: tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, float], row: int
    ) -> tuple[float, ...]:
        return (
            _binary64(item[0][block, row, field], label=f"{channel}.left"),
            _binary64(item[1][block, row, field], label=f"{channel}.left_rhs"),
            _binary64(item[2][block, row, field], label=f"{channel}.right"),
            _binary64(item[3][block, row, field], label=f"{channel}.right_rhs"),
            _binary64(item[4], label=f"{channel}.width"),
        )

    return tuple(
        (
            scalar(outer[0], row),
            tuple(scalar(item, row) for item in medium),
            tuple(scalar(item, row) for item in fine),
        )
        for row in range(OWNED_ROW_COUNT)
    )


def _row_line(ordinal: int, row: Sequence[object]) -> bytes:
    outer, medium, fine = row
    values: list[float] = []
    for segment in (outer, *medium, *fine):
        values.extend(_binary64(item, label=f"row[{ordinal}]") for item in segment)
    return (f"{ordinal}|" + "|".join(item.hex() for item in values) + "\n").encode(
        "ascii"
    )


def _feed_coefficients(hasher: object, ordinal: int, cubic: Sequence[Fraction]) -> None:
    line = (
        f"{ordinal}|"
        + "|".join(f"{item.numerator}/{item.denominator}" for item in cubic)
        + "\n"
    )
    hasher.update(line.encode("ascii"))


def _metadata(level: str, row: int, subinterval: int) -> dict[str, object]:
    return {"level": level, "row_index": row, "subinterval": subinterval}


def _survivor_key(item: object) -> tuple[object, ...]:
    return (
        item.polynomial_ordinal,
        item.location,
        item.location_ordinal,
        tuple(item.metadata),
    )


def _survivor_digest(items: Sequence[object]) -> str:
    answer = sha256(_SURVIVOR_KEY_HASH_DOMAIN)
    for ordinal, location, location_ordinal, metadata in sorted(
        _survivor_key(item) for item in items
    ):
        encoded = ",".join(f"{key}={value}" for key, value in metadata)
        answer.update(
            f"{ordinal}|{location}|{location_ordinal}|{encoded}\n".encode("ascii")
        )
    return answer.hexdigest()


def _stationary_digest(cubics: Sequence[object]) -> str:
    answer = sha256(_STATIONARY_COUNT_DOMAIN)
    for ordinal, cubic in enumerate(cubics):
        roots = _primary_stationary_intervals(
            cubic.coefficients, depth=PRIMARY_REFINEMENT_DEPTH
        )
        answer.update(f"{ordinal}|{len(roots)}\n".encode("ascii"))
    return answer.hexdigest()


def _overlap(a: Fraction, b: Fraction, c: Fraction, d: Fraction) -> bool:
    return max(a, c) <= min(b, d)


def _candidate_overlap(first: object, second: object) -> bool:
    return _overlap(
        first.parameter_lower,
        first.parameter_upper,
        second.parameter_lower,
        second.parameter_upper,
    ) and _overlap(
        first.absolute_lower,
        first.absolute_upper,
        second.absolute_lower,
        second.absolute_upper,
    )


def _localize_independently(
    coefficients: Sequence[
        tuple[tuple[Fraction, Fraction, Fraction, Fraction], dict[str, object]]
    ],
    *,
    level: str,
    coefficient_sha256: str,
    ceiling: int,
) -> _IndependentDifference:
    primary_cubics = tuple(
        LocalCubic(cubic, metadata) for cubic, metadata in coefficients
    )
    second_cubics = tuple(
        IndependentLocalCubicV2(cubic, metadata) for cubic, metadata in coefficients
    )
    try:
        primary = localize_absolute_maximum(
            primary_cubics,
            maximum_candidates=ceiling,
            refinement_depth=PRIMARY_REFINEMENT_DEPTH,
        )
        second = localize_absolute_maximum_independently_v2(
            second_cubics, maximum_candidates=ceiling
        )
    except RootIsolationInconclusive as exc:
        raise TDG10QA1PREF1Error(
            "PREF1_EXACT_RESOURCE_INCONCLUSIVE", f"{level}:{exc.reason}"
        ) from exc
    except RuntimeError as exc:
        raise TDG10QA1PREF1Error(
            "PREF1_EXACT_RESOURCE_INCONCLUSIVE", f"{level}:{exc}"
        ) from exc
    primary_stationary = _stationary_digest(primary_cubics)
    independent_stationary = second.stationary_count_stream_sha256
    first_items = {_survivor_key(item): item for item in primary.candidates}
    second_items = {_survivor_key(item): item for item in second.candidates}
    keys_equal = (
        len(first_items) == len(primary.candidates)
        and len(second_items) == len(second.candidates)
        and tuple(sorted(first_items)) == tuple(sorted(second_items))
    )
    route_facts = {
        "polynomial_count": primary.polynomial_count == second.polynomial_count,
        "candidate_count": primary.candidate_count == second.candidate_count,
        "classification": primary.classification == second.classification,
        "stationary_count": primary_stationary == independent_stationary,
        "survivor_keys": keys_equal,
        "survivor_intervals": keys_equal
        and all(
            _candidate_overlap(first_items[key], second_items[key])
            for key in first_items
        ),
        "global_interval": _overlap(
            primary.global_absolute_lower,
            primary.global_absolute_upper,
            second.global_absolute_lower,
            second.global_absolute_upper,
        ),
        "primary_tolerance_absent": primary.tolerance_used is False,
        "independent_tolerance_absent": second.tolerance_used is False,
        "independent_boundary_clipping_absent": second.boundary_clipping_used is False,
    }
    if not all(route_facts.values()):
        _stop("PREF1_EXACT_ROUTE_DISAGREEMENT", (level, route_facts))
    lower = max(primary.global_absolute_lower, second.global_absolute_lower)
    upper = min(primary.global_absolute_upper, second.global_absolute_upper)
    if lower < 0 or upper < lower:
        _stop("PREF1_EXACT_INTERVAL_DRIFT", level)
    return _IndependentDifference(
        lower=lower,
        upper=upper,
        polynomial_count=primary.polynomial_count,
        candidate_count=primary.candidate_count,
        co_maximizer_count=len(primary.candidates),
        localization_classification=primary.classification,
        coefficient_stream_sha256=coefficient_sha256,
        survivor_key_stream_sha256=_survivor_digest(primary.candidates),
        primary_stationary_count_stream_sha256=primary_stationary,
        independent_stationary_count_stream_sha256=independent_stationary,
        primary_evaluator_id=PRIMARY_LOCALIZER_ID,
        independent_evaluator_id=INDEPENDENT_LOCALIZER_ID,
        maximum_candidates=ceiling,
    )


def _assess_rows_independently(rows: Sequence[object]) -> _IndependentAdmission:
    row_hash = sha256(_ROW_HASH_DOMAIN)
    d01_hash = sha256(_COEFFICIENT_HASH_DOMAIN)
    d12_hash = sha256(_COEFFICIENT_HASH_DOMAIN)
    d01_coefficients: list[
        tuple[tuple[Fraction, Fraction, Fraction, Fraction], dict[str, object]]
    ] = []
    d12_coefficients: list[
        tuple[tuple[Fraction, Fraction, Fraction, Fraction], dict[str, object]]
    ] = []
    for row_index, row in enumerate(rows):
        outer, medium, fine = row
        if len(outer) != 5 or len(medium) != 2 or len(fine) != 4:
            _stop("PREF1_ROW_SCHEMA_DRIFT", row_index)
        row_hash.update(_row_line(row_index, row))
        outer_cubic = exact_hermite_coefficients(outer)
        medium_cubics = tuple(exact_hermite_coefficients(item) for item in medium)
        fine_cubics = tuple(exact_hermite_coefficients(item) for item in fine)
        outer_width = Fraction(
            *_binary64(outer[4], label="outer width").as_integer_ratio()
        )
        if any(
            Fraction(*_binary64(item[4], label="medium width").as_integer_ratio())
            != outer_width / 2
            for item in medium
        ) or any(
            Fraction(*_binary64(item[4], label="fine width").as_integer_ratio())
            != outer_width / 4
            for item in fine
        ):
            _stop("PREF1_ROW_WIDTH_DRIFT", row_index)
        for half, child in enumerate(medium_cubics):
            cubic = subtract_exact_cubics(
                restrict_exact_cubic_to_half(outer_cubic, half=half), child
            )
            _feed_coefficients(d01_hash, len(d01_coefficients), cubic)
            d01_coefficients.append((cubic, _metadata("D01", row_index, half)))
        for parent_index, parent in enumerate(medium_cubics):
            for half in (0, 1):
                child_index = 2 * parent_index + half
                cubic = subtract_exact_cubics(
                    restrict_exact_cubic_to_half(parent, half=half),
                    fine_cubics[child_index],
                )
                _feed_coefficients(d12_hash, len(d12_coefficients), cubic)
                d12_coefficients.append(
                    (cubic, _metadata("D12", row_index, child_index))
                )
    if len(rows) != OWNED_ROW_COUNT:
        _stop("PREF1_ROW_COUNT_DRIFT", len(rows))
    d01 = _localize_independently(
        d01_coefficients,
        level="D01",
        coefficient_sha256=d01_hash.hexdigest(),
        ceiling=PER_CHANNEL_MAXIMUM_CANDIDATES_D01,
    )
    d12 = _localize_independently(
        d12_coefficients,
        level="D12",
        coefficient_sha256=d12_hash.hexdigest(),
        ceiling=PER_CHANNEL_MAXIMUM_CANDIDATES_D12,
    )
    decision = classify_tdg6_channel(
        CertifiedMagnitudeInterval(d01.lower, d01.upper),
        CertifiedMagnitudeInterval(d12.lower, d12.upper),
    )
    left = TDG6_ORDER_SQUARED_MULTIPLIER * d12.upper**2
    right = d01.lower**2
    exact_zero = d01.upper == 0 and d12.upper == 0
    sufficient_pass = (not exact_zero) and left <= right
    sufficient_failure = (
        (not exact_zero)
        and not sufficient_pass
        and TDG6_ORDER_SQUARED_MULTIPLIER * d12.lower**2 > d01.upper**2
    )
    combined = sha256(
        _COMBINED_HASH_DOMAIN
        + f"{d01.coefficient_stream_sha256}\n{d12.coefficient_stream_sha256}\n".encode(
            "ascii"
        )
    ).hexdigest()
    return _IndependentAdmission(
        d01=d01,
        d12=d12,
        decision=decision,
        row_count=len(rows),
        row_stream_sha256=row_hash.hexdigest(),
        combined_coefficient_stream_sha256=combined,
        sufficient_pass_left=left,
        sufficient_pass_right=right,
        sufficient_contraction_pass=sufficient_pass,
        sufficient_contraction_failure=sufficient_failure,
        threshold_inconclusive=(
            not exact_zero and not sufficient_pass and not sufficient_failure
        ),
        maximum_candidates_D01=PER_CHANNEL_MAXIMUM_CANDIDATES_D01,
        maximum_candidates_D12=PER_CHANNEL_MAXIMUM_CANDIDATES_D12,
        refinement_depth=PRIMARY_REFINEMENT_DEPTH,
    )


def _assess_prepared_independently(prepared: object) -> _IndependentAssessment:
    width = _validate_prepared_shadow(prepared)
    channels = tuple(
        _IndependentChannel(
            channel=channel,
            evidence=_assess_rows_independently(
                _channel_rows(prepared, channel, width=width)
            ),
        )
        for channel in CHANNEL_ORDER
    )
    failed = tuple(
        item.channel for item in channels if not item.evidence.decision.admission_passed
    )
    return _IndependentAssessment(
        channel_evidence=channels,
        complete_admission_passed=not failed,
        failed_channels=failed,
    )


def _difference_payload(evidence: object, *, level: str) -> dict[str, Any]:
    maximum = (
        PER_CHANNEL_MAXIMUM_CANDIDATES_D01
        if level == "D01"
        else PER_CHANNEL_MAXIMUM_CANDIDATES_D12
    )
    polynomial_count = D01_POLYNOMIAL_COUNT if level == "D01" else D12_POLYNOMIAL_COUNT
    payload = {
        "candidate_count": _count(
            evidence.candidate_count, label=f"{level}.candidate_count"
        ),
        "co_maximizer_count": _count(
            evidence.co_maximizer_count, label=f"{level}.co_maximizer_count"
        ),
        "coefficient_stream_sha256": _digest(
            evidence.coefficient_stream_sha256,
            label=f"{level}.coefficient_stream_sha256",
        ),
        "independent_evaluator_id": evidence.independent_evaluator_id,
        "independent_stationary_count_stream_sha256": _digest(
            evidence.independent_stationary_count_stream_sha256,
            label=f"{level}.independent_stationary",
        ),
        "localization_classification": evidence.localization_classification,
        "lower": _encode_rational(evidence.lower),
        "maximum_candidates": evidence.maximum_candidates,
        "polynomial_count": evidence.polynomial_count,
        "primary_evaluator_id": evidence.primary_evaluator_id,
        "primary_stationary_count_stream_sha256": _digest(
            evidence.primary_stationary_count_stream_sha256,
            label=f"{level}.primary_stationary",
        ),
        "routes_agree": evidence.routes_agree,
        "survivor_key_stream_sha256": _digest(
            evidence.survivor_key_stream_sha256,
            label=f"{level}.survivor_key_stream_sha256",
        ),
        "upper": _encode_rational(evidence.upper),
    }
    if (
        payload["maximum_candidates"] != maximum
        or payload["polynomial_count"] != polynomial_count
        or payload["routes_agree"] is not True
        or payload["primary_evaluator_id"] != PRIMARY_EVALUATOR_ID
        or payload["independent_evaluator_id"] != INDEPENDENT_EVALUATOR_ID
        or payload["localization_classification"] not in _LOCALIZATION_CLASSES
        or payload["candidate_count"] > maximum
        or payload["co_maximizer_count"] > payload["candidate_count"]
    ):
        _stop("PREF1_DIFFERENCE_CONTRACT_DRIFT", level)
    return payload


def _reclassify_channel(
    *,
    channel: str,
    d01_lower: Fraction,
    d01_upper: Fraction,
    d12_lower: Fraction,
    d12_upper: Fraction,
) -> dict[str, Any]:
    try:
        outer = CertifiedMagnitudeInterval(d01_lower, d01_upper)
        finest = CertifiedMagnitudeInterval(d12_lower, d12_upper)
        decision = classify_tdg6_channel(outer, finest)
    except (TypeError, ValueError) as exc:
        raise TDG10QA1PREF1Error("PREF1_RECLASSIFY_DRIFT", channel) from exc
    left = TDG6_ORDER_SQUARED_MULTIPLIER * d12_upper**2
    right = d01_lower**2
    exact_zero = d01_upper == 0 and d12_upper == 0
    sufficient_holds = left <= right
    sufficient_pass = (not exact_zero) and sufficient_holds
    if d01_lower > 0:
        order_pass = three_halves_order_passes_squared(
            outer_lower_squared=d01_lower**2,
            finest_upper_squared=d12_upper**2,
        )
        if order_pass is not sufficient_holds:
            _stop("PREF1_SUFFICIENT_IDENTITY_DRIFT", channel)
    return {
        "decision": decision,
        "left": left,
        "right": right,
        "sufficient_condition_holds": sufficient_holds,
        "sufficient_contraction_pass": sufficient_pass,
    }


def _channel_replay(item: object) -> dict[str, Any]:
    channel = item.channel
    evidence = item.evidence
    decision = evidence.decision
    classified = _reclassify_channel(
        channel=channel,
        d01_lower=evidence.d01.lower,
        d01_upper=evidence.d01.upper,
        d12_lower=evidence.d12.lower,
        d12_upper=evidence.d12.upper,
    )
    independent = classified["decision"]
    if (
        independent != decision
        or classified["left"] != evidence.sufficient_pass_left
        or classified["right"] != evidence.sufficient_pass_right
        or classified["sufficient_contraction_pass"]
        is not evidence.sufficient_contraction_pass
        or evidence.row_count != OWNED_ROW_COUNT
        or evidence.maximum_candidates_D01 != PER_CHANNEL_MAXIMUM_CANDIDATES_D01
        or evidence.maximum_candidates_D12 != PER_CHANNEL_MAXIMUM_CANDIDATES_D12
        or evidence.refinement_depth != PRIMARY_REFINEMENT_DEPTH
        or evidence.evaluator_id != TDG10_EVALUATOR_ID
    ):
        _stop("PREF1_RECLASSIFY_DRIFT", channel)
    return {
        "D01": _difference_payload(evidence.d01, level="D01"),
        "D12": _difference_payload(evidence.d12, level="D12"),
        "admission_passed": decision.admission_passed,
        "candidate_evidence_available": True,
        "channel": channel,
        "classification": decision.classification,
        "combined_coefficient_stream_sha256": _digest(
            evidence.combined_coefficient_stream_sha256,
            label=f"{channel}.combined",
        ),
        "evaluator_id": TDG10_EVALUATOR_ID,
        "hash_evidence_available": True,
        "maximum_candidates_D01": PER_CHANNEL_MAXIMUM_CANDIDATES_D01,
        "maximum_candidates_D12": PER_CHANNEL_MAXIMUM_CANDIDATES_D12,
        "order_threshold_passed": decision.order_threshold_passed,
        "order_threshold_resolved": decision.order_threshold_resolved,
        "refinement_depth": PRIMARY_REFINEMENT_DEPTH,
        "route_evidence_available": True,
        "row_count": OWNED_ROW_COUNT,
        "row_stream_sha256": _digest(
            evidence.row_stream_sha256, label=f"{channel}.row_stream"
        ),
        "sufficient_condition": SUFFICIENT_CONDITION,
        "sufficient_condition_holds": classified["sufficient_condition_holds"],
        "sufficient_contraction_failure": evidence.sufficient_contraction_failure,
        "sufficient_contraction_pass": evidence.sufficient_contraction_pass,
        "sufficient_pass_left": _encode_rational(classified["left"]),
        "sufficient_pass_right": _encode_rational(classified["right"]),
        "temporal_retry_permitted": decision.temporal_retry_permitted,
        "threshold_inconclusive": evidence.threshold_inconclusive,
    }


def _compare_terminal_channel(live: Mapping[str, Any], raw: Mapping[str, Any]) -> None:
    channel = live["channel"]
    raw_upper = _fraction_from_rational(
        raw.get("d12_upper"), label=f"{channel}.raw.d12_upper"
    )
    live_upper = _fraction_from_rational(
        live["D12"]["upper"], label=f"{channel}.live.d12_upper"
    )
    if (
        raw.get("channel") != channel
        or raw.get("admission_passed") is not live["admission_passed"]
        or raw.get("classification") != live["classification"]
        or raw_upper != live_upper
    ):
        _stop("PREF1_TERMINAL_CHANNEL_DRIFT", channel)


def _recompute_channels(
    prepared: tdg6.TDG6PreparedGR0Compositor,
    raw_exact: Mapping[str, Any],
) -> dict[str, Any]:
    try:
        assessment = _assess_prepared_independently(prepared)
    except Exception as exc:
        if isinstance(exc, TDG10QA1PREF1Error):
            raise
        raise TDG10QA1PREF1Error("PREF1_EXACT_RECOMPUTE_DRIFT", exc) from exc
    raw_channels = raw_exact.get("channels")
    if not isinstance(raw_channels, list) or len(raw_channels) != len(CHANNEL_ORDER):
        _stop("PREF1_CHANNEL_COUNT_DRIFT", type(raw_channels).__name__)
    if tuple(item.channel for item in assessment.channel_evidence) != CHANNEL_ORDER:
        _stop("PREF1_CHANNEL_ORDER_DRIFT", "live channel order")
    channels: list[dict[str, Any]] = []
    failed: list[str] = []
    for item, raw_channel, (name, numerator, denominator) in zip(
        assessment.channel_evidence,
        raw_channels,
        EXPECTED_D12_UPPERS,
        strict=True,
    ):
        if item.channel != name:
            _stop("PREF1_CHANNEL_ORDER_DRIFT", (name, item.channel))
        if not isinstance(raw_channel, Mapping):
            _stop("PREF1_CHANNEL_TYPE_DRIFT", name)
        live = _channel_replay(item)
        _compare_terminal_channel(live, raw_channel)
        expected_upper = Fraction(int(numerator), int(denominator))
        observed_upper = _fraction_from_rational(
            live["D12"]["upper"], label=f"{name}.d12_upper"
        )
        if (
            observed_upper != expected_upper
            or live["classification"] != "resolved_order_pass"
            or live["admission_passed"] is not True
            or live["sufficient_contraction_pass"] is not True
            or live["sufficient_condition_holds"] is not True
        ):
            _stop("PREF1_CHANNEL_CLASS_DRIFT", name)
        if live["admission_passed"] is not True:
            failed.append(name)
        channels.append(live)
    complete = not failed
    if (
        assessment.complete_admission_passed is not complete
        or tuple(assessment.failed_channels) != tuple(failed)
        or raw_exact.get("complete_admission_passed") is not complete
        or tuple(raw_exact.get("failed_channels") or ()) != tuple(failed)
        or raw_exact.get("admission_is_all_of") is not True
        or raw_exact.get("interval_owner") != INTERVAL_OWNER
        or raw_exact.get("independent_route_required") is not True
        or raw_exact.get("channel_order") != list(CHANNEL_ORDER)
        or complete is not True
    ):
        _stop("PREF1_ALL_OF_DRIFT", failed)
    independent_class = reduce_qa1_terminal(
        complete_admission_passed=complete,
        failed_channels=failed,
    )
    return {
        **_expected_independent_replay_header(),
        "channel_summaries": _expected_channel_summaries(),
        "channels": channels,
        "raw_terminal_class_matches_independent_reduction": (
            independent_class == RAW_CLASSIFICATION
        ),
    }


def _decode_payload(raw: bytes) -> tuple[dict[str, np.ndarray], list[dict[str, Any]]]:
    if not isinstance(raw, bytes) or len(raw) > _MAX_RAW_LEAF_BYTES:
        _stop("PREF1_PAYLOAD_DRIFT", "archive byte budget")
    try:
        archive = zipfile.ZipFile(BytesIO(raw), "r")
        infos = archive.infolist()
    except (OSError, zipfile.BadZipFile) as exc:
        raise TDG10QA1PREF1Error("PREF1_PAYLOAD_DRIFT", "invalid ZIP") from exc
    expected_names = tuple(f"{name}.npy" for name in ARRAY_NAMES)
    try:
        observed_names = tuple(info.filename for info in infos)
        if observed_names != expected_names or len(set(observed_names)) != len(
            observed_names
        ):
            _stop("PREF1_PAYLOAD_DRIFT", observed_names)
        data: dict[str, np.ndarray] = {}
        manifest: list[dict[str, Any]] = []
        for name, info in zip(ARRAY_NAMES, infos, strict=True):
            mode = (info.external_attr >> 16) & 0o170000
            if (
                info.is_dir()
                or info.flag_bits & 0x1
                or info.file_size <= 0
                or info.file_size > _MAX_PAYLOAD_MEMBER_BYTES
                or Path(info.filename).name != info.filename
                or mode not in {0, stat.S_IFREG}
            ):
                _stop("PREF1_PAYLOAD_DRIFT", f"{name} archive entry is unsafe")
            encoded = archive.read(info)
            if len(encoded) != info.file_size:
                _stop("PREF1_PAYLOAD_DRIFT", f"{name} size differs")
            try:
                value = np.load(BytesIO(encoded), allow_pickle=False)
            except (OSError, ValueError) as exc:
                raise TDG10QA1PREF1Error("PREF1_PAYLOAD_DRIFT", name) from exc
            if (
                not isinstance(value, np.ndarray)
                or value.dtype != np.dtype("<f8")
                or not value.flags.c_contiguous
                or not value.shape
                or not np.isfinite(value).all()
            ):
                _stop("PREF1_PAYLOAD_DRIFT", f"{name} array differs")
            copied = value.copy(order="C")
            data[name] = copied
            digest = sha256(copied.tobytes(order="C")).hexdigest()
            manifest.append(
                {
                    "byte_count": int(copied.nbytes),
                    "bytes_sha256": digest,
                    "dtype": "<f8",
                    "name": name,
                    "order": "C",
                    "shape": list(copied.shape),
                }
            )
    finally:
        archive.close()
    if tuple(item["name"] for item in manifest) != ARRAY_NAMES:
        _stop("PREF1_PAYLOAD_DRIFT", "array order")
    if manifest != [dict(item) for item in EXPECTED_PAYLOAD_MANIFEST]:
        _stop("PREF1_PAYLOAD_MANIFEST_DRIFT", "per-array raw/semantic manifest")
    if sha256(_canonical(manifest)).hexdigest() != PAYLOAD_SEMANTIC_SHA256:
        _stop("PREF1_PAYLOAD_SEMANTIC_DRIFT", PAYLOAD_SEMANTIC_SHA256)
    return data, manifest


def _validate_descriptor(
    descriptor: Mapping[str, Any],
    *,
    descriptor_raw: bytes,
    payload_raw: bytes,
) -> dict[str, Any]:
    stored = descriptor.get("descriptor_sha256")
    body = {
        key: value for key, value in descriptor.items() if key != "descriptor_sha256"
    }
    content = sha256(_canonical(body)).hexdigest()
    if (
        stored != DESCRIPTOR_CONTENT_ADDRESS
        or content != DESCRIPTOR_CONTENT_ADDRESS
        or sha256(descriptor_raw).hexdigest() != DESCRIPTOR_RAW_SHA256
    ):
        _stop("PREF1_DESCRIPTOR_ADDRESS_DRIFT", content)
    if (
        descriptor.get("schema") != DIAGNOSTIC_ENDPOINT_SCHEMA
        or descriptor.get("classification") != DIAGNOSTIC_CLASSIFICATION
        or descriptor.get("artifact_id") != QA1_ARTIFACT_ID
        or descriptor.get("diagnostic_qualification_only") is not True
    ):
        _stop("PREF1_DESCRIPTOR_SCHEMA_DRIFT", descriptor.get("classification"))
    for name in ENDPOINT_FALSE_FLAGS:
        if descriptor.get(name) is not False:
            _stop("PREF1_ENDPOINT_FLAG_DRIFT", name)
    predecessor = descriptor.get("predecessor")
    if not isinstance(predecessor, Mapping) or predecessor.get("generation") != (
        PREDECESSOR_GENERATION
    ):
        _stop("PREF1_DESCRIPTOR_PREDECESSOR_DRIFT", predecessor)
    if predecessor.get("generation") == FORBIDDEN_RETRY3_GENERATION:
        _stop("PREF1_GENERATION_DRIFT", 10)
    if (
        descriptor.get("payload_raw_sha256") != PAYLOAD_RAW_SHA256
        or descriptor.get("payload_semantic_sha256") != PAYLOAD_SEMANTIC_SHA256
        or sha256(payload_raw).hexdigest() != PAYLOAD_RAW_SHA256
    ):
        _stop("PREF1_PAYLOAD_HASH_DRIFT", "descriptor/payload identity")
    arrays, manifest = _decode_payload(payload_raw)
    if descriptor.get("payload_manifest") != manifest:
        _stop("PREF1_PAYLOAD_MANIFEST_DRIFT", "descriptor manifest")
    if tuple(arrays) != ARRAY_NAMES:
        _stop("PREF1_PAYLOAD_DRIFT", tuple(arrays))
    return _expected_diagnostic_binding()


def _validate_terminal_exact(value: object) -> Mapping[str, Any]:
    if not isinstance(value, Mapping):
        _stop("PREF1_ADMISSION_TYPE_DRIFT", type(value).__name__)
    expected_keys = {
        "admission_is_all_of",
        "channel_count",
        "channel_order",
        "channels",
        "complete_admission_passed",
        "failed_channels",
        "independent_route_required",
        "interval_owner",
    }
    if set(value) != expected_keys:
        _stop("PREF1_SCHEMA_DRIFT", ("exact_complete_C", sorted(expected_keys)))
    return value


def bind_raw_result(config_raw: bytes, repository: Path) -> dict[str, Any]:
    expected = expected_bound_payload(config_raw)
    root = Path(os.path.abspath(os.fspath(repository)))
    metadata = root.lstat()
    if stat.S_ISLNK(metadata.st_mode) or not stat.S_ISDIR(metadata.st_mode):
        _stop("PREF1_REPOSITORY_UNSAFE", root)
    blobs, manifest, terminal, descriptor = _raw_snapshot(root)
    if manifest != _expected_manifest():
        _stop("PREF1_RAW_MANIFEST_DRIFT", "manifest schema differs")
    _authenticate_authority(root, manifest)
    before_store = _snapshot_store(root)
    sealed_store = (SEALED_STORE_LEAF_COUNT, SEALED_STORE_SHA256)
    if before_store != sealed_store:
        _stop("PREF1_STORE_IDENTITY_DRIFT", before_store)
    base = _expected_terminal_base()
    if set(terminal) != set(base) | _COMPLETED_EXTRA_KEYS:
        _stop("PREF1_TERMINAL_SCHEMA_DRIFT", sorted(terminal))
    if any(terminal.get(name) != value for name, value in base.items()):
        _stop("PREF1_TERMINAL_NONCLAIM_DRIFT", terminal.get("classification"))
    receipt = terminal.get("replay_receipt")
    if receipt != EXPECTED_REPLAY_RECEIPT:
        _stop("PREF1_REPLAY_RECEIPT_DRIFT", receipt)
    if (
        terminal.get("shadow_path_count") != 7
        or terminal.get("shadow_proposal_count") != 7
        or terminal.get("SSPRK3_stage_and_endpoint_record_count") != 28
    ):
        _stop("PREF1_SHADOW_COUNT_DRIFT", "completed terminal counts")
    raw_exact = _validate_terminal_exact(terminal.get("exact_complete_C"))
    prepared, live_receipt = _restore_shadow(root)
    replay = _recompute_channels(prepared, raw_exact)
    independent_class = reduce_qa1_terminal(
        complete_admission_passed=True,
        failed_channels=(),
    )
    if (
        independent_class != RAW_CLASSIFICATION
        or terminal.get("classification") != independent_class
        or live_receipt != EXPECTED_REPLAY_RECEIPT
    ):
        _stop("PREF1_TERMINAL_CLASS_DRIFT", independent_class)
    diagnostic = _validate_descriptor(
        descriptor,
        descriptor_raw=blobs[DESCRIPTOR_RELATIVE],
        payload_raw=blobs[PAYLOAD_RELATIVE],
    )
    endpoint = terminal.get("diagnostic_fine_endpoint")
    if not isinstance(endpoint, Mapping):
        _stop("PREF1_ENDPOINT_TYPE_DRIFT", type(endpoint).__name__)
    if (
        endpoint.get("present") is not True
        or endpoint.get("accepted_state") is not False
        or endpoint.get("descriptor_relative") != DESCRIPTOR_RELATIVE
        or endpoint.get("descriptor_sha256") != DESCRIPTOR_CONTENT_ADDRESS
        or endpoint.get("payload_relative") != PAYLOAD_RELATIVE
        or endpoint.get("payload_raw_sha256") != PAYLOAD_RAW_SHA256
        or endpoint.get("payload_semantic_sha256") != PAYLOAD_SEMANTIC_SHA256
        or endpoint.get("payload_manifest")
        != [dict(item) for item in EXPECTED_PAYLOAD_MANIFEST]
    ):
        _stop("PREF1_ENDPOINT_RECORD_DRIFT", "terminal diagnostic record")
    after_store = _snapshot_store(root)
    blobs_after = _raw_tree(root)
    if after_store != before_store or blobs_after != blobs:
        _stop("PREF1_INPUT_MUTATED_DURING_BIND", "raw/store input changed")
    return {
        **expected,
        "diagnostic_binding": diagnostic,
        "independent_replay": replay,
    }


def build_pref1_result(
    config_raw: bytes, repository: Path, *, live: bool = True
) -> dict[str, Any]:
    payload = (
        bind_raw_result(config_raw, repository)
        if live
        else expected_evidence(config_raw)
    )
    return {"artifact_id": ARTIFACT_ID, "artifact_payload": payload}


def _validate_difference(
    value: object, *, level: str, expected_upper: Fraction | None = None
) -> tuple[Fraction, Fraction]:
    payload = _mapping(value, _DIFFERENCE_KEYS, f"{level}")
    lower = _fraction_from_rational(payload.get("lower"), label=f"{level}.lower")
    upper = _fraction_from_rational(payload.get("upper"), label=f"{level}.upper")
    if lower < 0 or upper < 0 or lower > upper:
        _stop("PREF1_INTERVAL_DRIFT", level)
    maximum = (
        PER_CHANNEL_MAXIMUM_CANDIDATES_D01
        if level.endswith("D01") or level == "D01"
        else PER_CHANNEL_MAXIMUM_CANDIDATES_D12
    )
    polynomial_count = (
        D01_POLYNOMIAL_COUNT
        if level.endswith("D01") or level == "D01"
        else D12_POLYNOMIAL_COUNT
    )
    candidate_count = _count(
        payload.get("candidate_count"), label=f"{level}.candidate_count"
    )
    co_maximizers = _count(
        payload.get("co_maximizer_count"), label=f"{level}.co_maximizer_count"
    )
    if (
        payload.get("maximum_candidates") != maximum
        or payload.get("polynomial_count") != polynomial_count
        or payload.get("routes_agree") is not True
        or payload.get("primary_evaluator_id") != PRIMARY_EVALUATOR_ID
        or payload.get("independent_evaluator_id") != INDEPENDENT_EVALUATOR_ID
        or payload.get("localization_classification") not in _LOCALIZATION_CLASSES
        or candidate_count > maximum
        or co_maximizers > candidate_count
    ):
        _stop("PREF1_DIFFERENCE_CONTRACT_DRIFT", level)
    for name in (
        "coefficient_stream_sha256",
        "independent_stationary_count_stream_sha256",
        "primary_stationary_count_stream_sha256",
        "survivor_key_stream_sha256",
    ):
        _digest(payload.get(name), label=f"{level}.{name}")
    if expected_upper is not None and upper != expected_upper:
        _stop("PREF1_D12_UPPER_DRIFT", level)
    return lower, upper


def _validate_rich_channel(
    value: object, *, expected: tuple[str, str, str]
) -> dict[str, Any]:
    channel, numerator, denominator = expected
    payload = _mapping(value, _CHANNEL_REPLAY_KEYS, channel)
    if payload.get("channel") != channel:
        _stop("PREF1_CHANNEL_ORDER_DRIFT", (channel, payload.get("channel")))
    expected_upper = Fraction(int(numerator), int(denominator))
    d01_lower, d01_upper = _validate_difference(
        payload.get("D01"), level=f"{channel}.D01"
    )
    d12_lower, d12_upper = _validate_difference(
        payload.get("D12"), level=f"{channel}.D12", expected_upper=expected_upper
    )
    classified = _reclassify_channel(
        channel=channel,
        d01_lower=d01_lower,
        d01_upper=d01_upper,
        d12_lower=d12_lower,
        d12_upper=d12_upper,
    )
    decision = classified["decision"]
    left = _fraction_from_rational(
        payload.get("sufficient_pass_left"), label=f"{channel}.left"
    )
    right = _fraction_from_rational(
        payload.get("sufficient_pass_right"), label=f"{channel}.right"
    )
    if (
        decision.classification != payload.get("classification")
        or decision.admission_passed is not payload.get("admission_passed")
        or decision.temporal_retry_permitted
        is not payload.get("temporal_retry_permitted")
        or decision.order_threshold_resolved
        is not payload.get("order_threshold_resolved")
        or decision.order_threshold_passed is not payload.get("order_threshold_passed")
        or classified["left"] != left
        or classified["right"] != right
        or classified["sufficient_condition_holds"]
        is not payload.get("sufficient_condition_holds")
        or classified["sufficient_contraction_pass"]
        is not payload.get("sufficient_contraction_pass")
        or payload.get("classification") != "resolved_order_pass"
        or payload.get("admission_passed") is not True
        or payload.get("sufficient_contraction_pass") is not True
        or payload.get("sufficient_condition") != SUFFICIENT_CONDITION
        or payload.get("row_count") != OWNED_ROW_COUNT
        or payload.get("maximum_candidates_D01") != PER_CHANNEL_MAXIMUM_CANDIDATES_D01
        or payload.get("maximum_candidates_D12") != PER_CHANNEL_MAXIMUM_CANDIDATES_D12
        or payload.get("refinement_depth") != PRIMARY_REFINEMENT_DEPTH
        or payload.get("evaluator_id") != TDG10_EVALUATOR_ID
        or payload.get("route_evidence_available") is not True
        or payload.get("candidate_evidence_available") is not True
        or payload.get("hash_evidence_available") is not True
    ):
        _stop("PREF1_RECLASSIFY_DRIFT", channel)
    _digest(payload.get("row_stream_sha256"), label=f"{channel}.row_stream")
    _digest(
        payload.get("combined_coefficient_stream_sha256"),
        label=f"{channel}.combined",
    )
    return payload


def _validate_independent_replay(value: object) -> None:
    if not isinstance(value, Mapping):
        _stop("PREF1_COMPACT_DRIFT", "independent replay absent")
    header = _expected_independent_replay_header()
    if set(value) != set(header) | {"channel_summaries", "channels"}:
        _stop("PREF1_COMPACT_DRIFT", "independent replay schema")
    for key, expected in header.items():
        if value.get(key) != expected:
            _stop("PREF1_COMPACT_DRIFT", f"independent replay {key}")
    summaries = value.get("channel_summaries")
    if summaries != _expected_channel_summaries():
        _stop("PREF1_COMPACT_DRIFT", "channel summaries")
    channels = value.get("channels")
    if not isinstance(channels, list) or len(channels) != len(CHANNEL_ORDER):
        _stop("PREF1_COMPACT_DRIFT", "rich channels absent")
    failed: list[str] = []
    for item, expected in zip(channels, EXPECTED_D12_UPPERS, strict=True):
        live = _validate_rich_channel(item, expected=expected)
        if live["admission_passed"] is not True:
            failed.append(expected[0])
    if failed or value.get("complete_admission_passed") is not True:
        _stop("PREF1_ALL_OF_DRIFT", failed)
    independent_class = reduce_qa1_terminal(
        complete_admission_passed=True,
        failed_channels=(),
    )
    if (
        independent_class != RAW_CLASSIFICATION
        or value.get("raw_terminal_class_matches_independent_reduction") is not True
    ):
        _stop("PREF1_TERMINAL_CLASS_DRIFT", independent_class)


def validate_compact_result(config_raw: bytes, result_raw: bytes) -> dict[str, Any]:
    expected = expected_compact_result(config_raw)
    if sha256(result_raw).hexdigest() != COMPACT_RESULT_SHA256:
        _stop("PREF1_COMPACT_HASH_DRIFT", COMPACT_RESULT_SHA256)
    result = _json(result_raw, "compact result")
    if (
        set(result) != {"artifact_id", "artifact_payload"}
        or result.get("artifact_id") != ARTIFACT_ID
        or not isinstance(result.get("artifact_payload"), Mapping)
    ):
        _stop("PREF1_COMPACT_DRIFT", "compact identity differs")
    payload = result["artifact_payload"]
    expected_payload = expected["artifact_payload"]
    if set(payload) != set(expected_payload):
        _stop("PREF1_COMPACT_DRIFT", "compact payload schema")
    for key, value in expected_payload.items():
        if key == "independent_replay":
            continue
        if payload.get(key) != value:
            _stop("PREF1_COMPACT_DRIFT", key)
    _validate_independent_replay(payload.get("independent_replay"))
    return result


__all__ = (
    "ARTIFACT_ID",
    "COMPACT_RESULT_SHA256",
    "CONFIG_PATH",
    "OWNER_DOCUMENT",
    "RESULT_PATH",
    "TDG10QA1PREF1Error",
    "bind_raw_result",
    "build_pref1_result",
    "canonical_result",
    "expected_bound_payload",
    "expected_compact_result",
    "expected_config",
    "expected_evidence",
    "reduce_qa1_terminal",
    "validate_compact_result",
)
