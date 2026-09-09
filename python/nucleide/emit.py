"""Single-material emission to legacy transport-code cards (backed by the `nucleide-emit` crate)."""

from nucleide._internal import (
    emit_armi_cards,
    emit_armi_drift_table,
    emit_cards,
    emit_drift_table,
)

__all__ = [
    "emit_armi_cards",
    "emit_armi_drift_table",
    "emit_cards",
    "emit_drift_table",
]
