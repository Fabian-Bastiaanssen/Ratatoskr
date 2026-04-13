import shutil 

from loguru import logger
from tqdm import tqdm
from pathlib import Path
from ratatoskr.utils import get_credentials, make_dir, get_version
from ratatoskr.lpsn import set_lpsn_client
from ratatoskr.bacdive import set_bacdive_client

def set_logger_format(debug: bool = False) -> str:
    """
    Set the logger format for console output with color-coded timestamps and levels.
    
    Removes the default logger format and configures a custom format that outputs
    to tqdm-compatible console with magenta timestamps, colored log levels, and
    black message text.
    
    Returns:
        str: The formatted logger string pattern used for logging.
    """
    
    logger_format = "<magenta>{time:YYYY-MM-DD HH:mm:ss}</magenta> | <level>{level: <8}</level> | <level>{message}</level>"
    logger.remove() # remove default logger format
    if debug:
        logger.add(lambda msg: tqdm.write(msg, end=""), colorize=True, format= logger_format, level="DEBUG")
    else:
        logger.add(lambda msg: tqdm.write(msg, end=""), colorize=True, format= logger_format, level="INFO")

    return logger_format


def create_log_file(output_path: Path, logger_format: str) -> None:
    """
    Create a log file in the specified output directory.
    
    Configures file-based logging with rotation and compression settings.
    The log file is created at output_path/ratatoskr.log with automatic rotation
    at 100 MB and zip compression.
    
    Args:
        output_path: The directory path where the log file will be created.
        logger_format: The format string to use for log entries.
    """
    logger.add(
        output_path / "ratatoskr.log",
        format=logger_format,
        level=1,
        rotation="100 MB",
        compression="zip",
    )

def validate_out_folder(output_path: Path, force: bool) -> None:
    """
    Validate and prepare the output folder.
    
    Checks if the output directory exists. If it does and force is True, removes
    the existing directory. If it exists and force is False, logs an error and
    exits. Finally, creates the output directory.
    
    Args:
        output_path: The path to the output directory to validate.
        force: Whether to forcefully overwrite an existing output directory.
        
    Raises:
        SystemExit: If output directory exists and force is False.
    """

    if output_path.exists():
        if force is True:
            shutil.rmtree(output_path)
        else:
            logger.error(
                    "Output dir exists and --force not specified; ❌ failed check." 
                )
            logger.info(
                    "Please specify -f or --force to overwrite the output directory, or alternatively specify a different output directory."
                )
            sys.exit(1)

    make_dir(output_path)


def initialise_clients(dev_mode=False):

    email, password, _ = get_credentials(email=True, password=True, api_key=False, api_being_accessed="LPSN", dev_mode=dev_mode)

    lpsn_client = set_lpsn_client(email, password)
    bacdive_client = set_bacdive_client(email, password)

    return lpsn_client, bacdive_client


def set_up_logger(output_path: Path, force: bool, debug: bool = False) -> None:
    """
    Set up the complete logging system for ratatoskr.
    
    Initializes the logger with custom formatting, validates the output folder,
    creates the log file, and logs the ratatoskr version header.
    
    Args:
        output_path: The directory path where logs will be stored.
        force: Whether to forcefully overwrite an existing output directory.
    """
    
    logger_format = set_logger_format(debug=debug)
    validate_out_folder(output_path, force)
    create_log_file(output_path, logger_format)
    
    logger.info("###################################")
    logger.info("#####    Ratatoskr v" + get_version() + "     #####")
    logger.info("###################################\n")
    # provide current date, time and version info in log header for better tracking and reproducibility
    logger.info(f"Ran on {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}, with ratatoskr version {get_version()}")
