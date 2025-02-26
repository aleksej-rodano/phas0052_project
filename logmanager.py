import logging
import os


class LogManager:
    def _init_(self, log_file="system.log"):
        """
        Inizializza il LogManager con un file di log.
        :param log_file: Nome del file di log (default: 'system.log').
        """
        self.log_file = log_file
        self.setup_logger()

    def setup_logger(self):
        """
        Configura il logger.
        """
        self.logger = logging.getLogger("SystemLogger")
        self.logger.setLevel(logging.DEBUG)

        # Formattazione del log
        formatter = logging.Formatter("%(asctime)s - %(levelname)s - %(message)s")

        # File Handler per scrivere i log su file
        file_handler = logging.FileHandler(self.log_file)
        file_handler.setLevel(logging.DEBUG)
        file_handler.setFormatter(formatter)

        # Stream Handler (opzionale) per stampare i log anche nel terminale
        console_handler = logging.StreamHandler()
        console_handler.setLevel(logging.INFO)
        console_handler.setFormatter(formatter)

        # Aggiungi gli handler al logger
        if not self.logger.hasHandlers():
            self.logger.addHandler(file_handler)
            self.logger.addHandler(console_handler)

    def log(self, message, level="info"):
        """
        Registra un messaggio nel log.
        :param message: Messaggio da loggare.
        :param level: Livello del log ('info', 'debug', 'warning', 'error', 'critical').
        """
        if level == "info":
            self.logger.info(message)
        elif level == "debug":
            self.logger.debug(message)
        elif level == "warning":
            self.logger.warning(message)
        elif level == "error":
            self.logger.error(message)
        elif level == "critical":
            self.logger.critical(message)
        else:
            self.logger.info(message)

    def clear_log(self):
        """
        Cancella il contenuto del file di log.
        """
        try:
            with open(self.log_file, "w") as file:
                file.write("")
            self.log("Log file cleared.", "info")
        except Exception as e:
            self.log(f"Error clearing log file: {e}", "error")


if _name_ == "_main_":
    # Esempio d'uso
    log_manager = LogManager("system.log")
