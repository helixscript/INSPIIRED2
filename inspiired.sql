CREATE TABLE fragments (
    trial VARCHAR(100) NOT NULL,
    subject VARCHAR(100) NOT NULL,
    sample VARCHAR(100) NOT NULL,
    replicate INT NOT NULL,
    ref_genome VARCHAR(10) NOT NULL,
    mode VARCHAR(10) NOT NULL,
    total_fragments INT,
    processed_date TIMESTAMP DEFAULT CURRENT_TIMESTAMP,
    data_file_name VARCHAR(50),
    PRIMARY KEY (trial, subject, sample, replicate, ref_genome, mode)
);