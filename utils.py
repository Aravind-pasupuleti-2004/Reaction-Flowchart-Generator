import re
import logging
from datetime import datetime

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(levelname)s - %(message)s",
    handlers=[
        logging.FileHandler("conversion.log", encoding="utf-8"),
        logging.StreamHandler(),
    ],
)
logger = logging.getLogger(__name__)


def detect_input_type(text: str) -> str:
    if not text or not text.strip():
        return "unknown"

    text = text.strip()

    if re.match(r"^InChI=1S?/", text):
        return "inchi"

    if re.match(r"^InChIKey=", text):
        return "inchikey"

    if re.match(r"^[A-Za-z0-9@+\-\[\]()\\\/%#$.=,:;]+$", text):
        has_organic = bool(re.search(r"[BCNOPSFIHbru]", text))
        has_branches = "(" in text or ")" in text
        has_rings = "1" in text or "2" in text
        if (has_organic or has_branches or has_rings) and len(text) > 1:
            try:
                from rdkit import RDLogger, Chem
                rd_logger = RDLogger.logger()
                rd_logger.setLevel(RDLogger.FATAL)
                mol = Chem.MolFromSmiles(text)
                rd_logger.setLevel(RDLogger.WARNING)
                if mol is not None:
                    return "smiles"
            except Exception:
                pass

    if " " in text or text[0].isupper():
        return "iupac_or_common"

    return "iupac_or_common"


def detect_file_type(file_path: str) -> str:
    if file_path.endswith(".xlsx"):
        return "excel"
    elif file_path.endswith(".csv"):
        return "csv"
    return "unknown"


def sanitize_filename(name: str) -> str:
    return re.sub(r'[\\/*?:"<>|]', "_", name)


def timestamp() -> str:
    return datetime.now().strftime("%Y%m%d_%H%M%S")
