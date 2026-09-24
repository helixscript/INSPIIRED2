CREATE DATABASE IF NOT EXISTS inspiired2;
USE inspiired2;

CREATE TABLE fragments (
    trial VARCHAR(100) NOT NULL,
    subject VARCHAR(100) NOT NULL,
    sample VARCHAR(100) NOT NULL,
    replicate INT NOT NULL,
    ref_genome VARCHAR(10) NOT NULL,
    mode VARCHAR(20) NOT NULL,
    total_fragments INT,
    processed_date TIMESTAMP DEFAULT CURRENT_TIMESTAMP,
    data_file_name VARCHAR(50),
    PRIMARY KEY (trial, subject, sample, replicate, ref_genome, mode)
) ENGINE=InnoDB;

CREATE TABLE sites (
    trial VARCHAR(100) NOT NULL,
    subject VARCHAR(100) NOT NULL,
    sample VARCHAR(100) NOT NULL,
    ref_genome VARCHAR(10) NOT NULL,
    mode VARCHAR(20) NOT NULL,
    total_sites INT,
    processed_date TIMESTAMP DEFAULT CURRENT_TIMESTAMP,
    data_file_name VARCHAR(50),
    PRIMARY KEY (trial, subject, sample, ref_genome, mode)
) ENGINE=InnoDB;

-- Permanent faux records allow read-only access tests before real data is loaded.
INSERT INTO fragments
    (trial, subject, sample, replicate, ref_genome, mode,
     total_fragments, processed_date, data_file_name)
VALUES
    ('fauxtrial', 'fauxSubject', 'fauxSample', 0, 'xxx', 'xxx',
     0, '2000-01-01 00:00:00', 'xxx');

INSERT INTO sites
    (trial, subject, sample, ref_genome, mode,
     total_sites, processed_date, data_file_name)
VALUES
    ('fauxtrial', 'fauxSubject', 'fauxSample', 'xxx', 'xxx',
     0, '2000-01-01 00:00:00', 'xxx');

-- Create accounts and assign database privileges.
CREATE USER 'inspiired2_user'@'%' IDENTIFIED BY 'user@+1';
CREATE USER 'inspiired2_admin'@'%' IDENTIFIED BY 'admin@+2';

GRANT SELECT ON `inspiired2`.*
TO 'inspiired2_user'@'%';

GRANT SELECT, INSERT, UPDATE, DELETE ON `inspiired2`.*
TO 'inspiired2_admin'@'%';
